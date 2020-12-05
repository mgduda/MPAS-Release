! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module test_core_timekeeping_tests

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_field_routines
   use mpas_dmpar
   use mpas_threading
   use mpas_log
   use mpas_timer
   use mpas_timekeeping

   implicit none
   private

   public :: test_core_test_intervals, &
             test_core_test_alarms

   contains

   !***********************************************************************
   !
   !  routine test_core_test_intervals
   !
   !> \brief   MPAS Test Core timekeeping interval tests
   !> \author  Doug Jacobsen
   !> \date    08/23/2016
   !> \details 
   !>  This routine performs tests related to timekeeping intervals
   !
   !-----------------------------------------------------------------------
   subroutine test_core_test_intervals(domain, threadErrs, err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      character (len=StrKIND) :: int1_str, int2_str
      character (len=StrKIND) :: ref_str

      integer :: threadNum
      integer :: iErr, err_tmp
      integer :: divs

      call mpas_timer_start('timekeeping tests')

      threadNum = mpas_threading_get_thread_num()
      err = 0

      if ( threadNum == 0 ) then
         call mpas_log_write(' Performing time interval tests')

         call mpas_log_write('  Test 1:')
         call test_core_interval_test('0001-01-01_00:00:00', '0000-01-00_10:00:00', '0001_00:00:00', 31_I8KIND, '0000_10:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            call mpas_log_write('   Result: PASSED')
         else
            call mpas_log_write(' * Result: FAILED', MPAS_LOG_ERR)
         end if

         call mpas_log_write('  Test 2:')
         call test_core_interval_test('0001-01-01_00:00:00', '0000-01-00_00:00:00', '0001_00:00:00', 31_I8KIND, '0000_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            call mpas_log_write('   Result: PASSED')
         else
            call mpas_log_write(' * Result: FAILED', MPAS_LOG_ERR)
         end if

         call mpas_log_write('  Test 3:')
         call test_core_interval_test('0001-02-01_00:00:00', '0000-01-00_10:00:00', '0001_00:00:00', 28_I8KIND, '0000_10:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            call mpas_log_write('   Result: PASSED')
         else
            call mpas_log_write(' * Result: FAILED', MPAS_LOG_ERR)
         end if

         call mpas_log_write('  Test 4:')
         call test_core_interval_test('0001-02-01_00:00:00', '0000-01-00_00:00:00', '0001_00:00:00', 28_I8KIND, '0000_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            call mpas_log_write('   Result: PASSED')
         else
            call mpas_log_write(' * Result: FAILED', MPAS_LOG_ERR)
         end if

         call mpas_log_write('  Test 5:')
         call test_core_interval_test('0001-01-01_00:00:00', '0000-00-00_01:00:00', '0000_00:30:00', 2_I8KIND, '0000_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            call mpas_log_write('   Result: PASSED')
         else
            call mpas_log_write(' * Result: FAILED', MPAS_LOG_ERR)
         end if

         call mpas_log_write('  Test 6:')
         call test_core_interval_test('0001-01-01_00:00:00', '0001-01-00_00:00:00', '0001-00-00_00:00:00', 1_I8KIND, '0000-00-31_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            call mpas_log_write('   Result: PASSED')
         else
            call mpas_log_write(' * Result: FAILED', MPAS_LOG_ERR)
         end if

         call mpas_log_write('  Test 7:')
         call test_core_interval_test('0000-01-01_00:00:00', '1850-00-00_00:00:00', '00:00:01', 58341600000_I8KIND, '0000-00-00_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            call mpas_log_write('   Result: PASSED')
         else
            call mpas_log_write(' * Result: FAILED', MPAS_LOG_ERR)
         end if

         call mpas_log_write(' Completed time interval tests')

      end if

      call mpas_timer_stop('timekeeping tests')

   end subroutine test_core_test_intervals!}}}


   !***********************************************************************
   !
   !  routine test_core_test_alarms
   !
   !> \brief   MPAS Test Core timekeeping alarm tests
   !> \author  Michael Duda
   !> \date    3 December 2020
   !> \details 
   !>  Cases to test:
   !>   1) Alarms are added and removed
   !>   2) The current time on the clock is changed
   !>   3) Alarms before the current time are ignored
   !>   4) The clock direction is changed
   !
   !-----------------------------------------------------------------------
   subroutine test_core_test_alarms(domain, threadErrs, err)

      use mpas_derived_types, only : MPAS_Clock_type, MPAS_TimeInterval_type, MPAS_Time_type
      use mpas_timekeeping, only : mpas_create_clock, mpas_destroy_clock, mpas_set_time, mpas_set_timeInterval

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      integer :: i
      integer :: threadNum
      integer :: ierr
      type (MPAS_Clock_type) :: clock
      type (MPAS_Time_type) :: now
      type (MPAS_Time_type) :: startTime
      type (MPAS_Time_type) :: stopTime
      type (MPAS_Time_type) :: alarmTime
      type (MPAS_Time_type) :: nextRingTime
      type (MPAS_TimeInterval_type) :: timeStep
      type (MPAS_TimeInterval_type) :: alarmInterval
      character(len=StrKIND) :: temp

      call mpas_timer_start('alarm tests')

      threadNum = mpas_threading_get_thread_num()
      err = 0
      threadErrs(:) = 0

      if (threadNum == 0) then
         call mpas_set_time(startTime, dateTimeString='1850-01-01_00:00:00', ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error setting startTime in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if

         call mpas_set_time(stopTime, dateTimeString='2000-01-01_00:00:00', ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error setting startTime in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if

         call mpas_set_timeInterval(timeStep, timeString='01:00:00', ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error setting timeStep in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if

         call mpas_create_clock(clock, startTime, timeStep, stopTime=stopTime, ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error creating clock in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if


         !
         ! Add a recurring alarm to the clock
         !
         call mpas_set_time(alarmTime, dateTimeString='1850-01-01_00:00:00', ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error setting alarmTime in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if

         call mpas_set_timeInterval(alarmInterval, timeString='24:00:00', ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error setting timeStep in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if

         call mpas_add_clock_alarm(clock, 'daily', alarmTime, alarmTimeInterval=alarmInterval, ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error adding daily alarm in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if


         !
         ! Add a one-shot alarm to the clock
         !
         call mpas_set_time(alarmTime, dateTimeString='1850-01-02_00:00:00', ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error setting alarmTime in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if

         call mpas_add_clock_alarm(clock, 'one-shot', alarmTime, ierr=ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error adding one-shot alarm in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if

         nextRingTime = mpas_alarm_next_ring_time(clock, 'daily', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time: '//trim(temp))

         nextRingTime = mpas_alarm_next_ring_time(clock, 'one-shot', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time: '//trim(temp))

         if (mpas_is_alarm_ringing(clock, 'daily')) then
            call mpas_reset_clock_alarm(clock, 'daily')
         end if

         if (mpas_is_alarm_ringing(clock, 'one-shot')) then
            call mpas_reset_clock_alarm(clock, 'one-shot')
         end if

         !
         ! Directly set the time 12 hours ahead
         !
         call mpas_set_time(now, dateTimeString='1850-01-01_12:00:00')
         call mpas_set_clock_time(clock, now, MPAS_NOW)

         now = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(now, dateTimeString=temp)
         call mpas_log_write('Current time: '//trim(temp))

         nextRingTime = mpas_alarm_next_ring_time(clock, 'daily', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time (daily): '//trim(temp))

         nextRingTime = mpas_alarm_next_ring_time(clock, 'one-shot', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time (one-shot): '//trim(temp))

         !
         ! Change the direction of the clock
         !
         call mpas_log_write('Clock is now running backwards')
         call mpas_set_clock_direction(clock, MPAS_BACKWARD)

         nextRingTime = mpas_alarm_next_ring_time(clock, 'daily', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time (daily): '//trim(temp))

         nextRingTime = mpas_alarm_next_ring_time(clock, 'one-shot', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time (one-shot): '//trim(temp))

         !
         ! Directly set the time 36 hours ahead
         !
         call mpas_set_time(now, dateTimeString='1850-01-02_12:00:00')
         call mpas_set_clock_time(clock, now, MPAS_NOW)

         now = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(now, dateTimeString=temp)
         call mpas_log_write('Current time: '//trim(temp))

         nextRingTime = mpas_alarm_next_ring_time(clock, 'daily', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time (daily): '//trim(temp))

         nextRingTime = mpas_alarm_next_ring_time(clock, 'one-shot', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time (one-shot): '//trim(temp))

         !
         ! Change the direction of the clock
         !
         call mpas_log_write('Clock is now running forwards')
         call mpas_set_clock_direction(clock, MPAS_FORWARD)

         nextRingTime = mpas_alarm_next_ring_time(clock, 'daily', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time (daily): '//trim(temp))

         nextRingTime = mpas_alarm_next_ring_time(clock, 'one-shot', strictlyInFuture=.true., ierr=ierr)
         call mpas_get_time(nextRingTime, dateTimeString=temp)
         call mpas_log_write('Next ring time (one-shot): '//trim(temp))

         !
         ! Directly set the time back to start time
         !
         now = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
         call mpas_set_clock_time(clock, now, MPAS_NOW)

         !
         ! Advance clock by 48 time steps = two days
         !
         do i = 1, 48

             call mpas_advance_clock(clock, ierr=ierr)
             if (ierr /= 0) then
                 call mpas_log_write('Error advancing clock in test_core_test_alarms', messageType=MPAS_LOG_ERR)
                 err = 1
                 threadErrs(1) = 1
                 return
             end if

             now = mpas_get_clock_time(clock, MPAS_NOW, ierr)
             call mpas_get_time(now, dateTimeString=temp)
             call mpas_log_write('Current time: '//trim(temp))

             nextRingTime = mpas_alarm_next_ring_time(clock, 'daily', strictlyInFuture=.true., ierr=ierr)
             call mpas_get_time(nextRingTime, dateTimeString=temp)
             call mpas_log_write('Next ring time (daily): '//trim(temp))

             nextRingTime = mpas_alarm_next_ring_time(clock, 'one-shot', strictlyInFuture=.true., ierr=ierr)
             call mpas_get_time(nextRingTime, dateTimeString=temp)
             call mpas_log_write('Next ring time (one-shot): '//trim(temp))

             call mpas_log_write('-------------')

             if (mpas_is_alarm_ringing(clock, 'daily')) then
                call mpas_reset_clock_alarm(clock, 'daily')
             end if

             if (mpas_is_alarm_ringing(clock, 'one-shot')) then
                call mpas_reset_clock_alarm(clock, 'one-shot')
             end if
         end do

         call mpas_destroy_clock(clock, ierr)
         if (ierr /= 0) then
             call mpas_log_write('Error destroying clock in test_core_test_alarms', messageType=MPAS_LOG_ERR)
             err = 1
             threadErrs(1) = 1
             return
         end if

         call mpas_log_write(' Completed alarm tests')
      end if

      call mpas_timer_stop('alarm tests')

   end subroutine test_core_test_alarms


   subroutine test_core_interval_test(ref_str, int1_str, int2_str, expected_divs, expected_remainder_str, ierr)!{{{
      character (len=*), intent(in) :: ref_str, int1_str, int2_str
      integer (kind=I8KIND), intent(in) :: expected_divs
      character (len=*), intent(in) :: expected_remainder_str
      integer, intent(out) :: ierr

      integer (kind=I8KIND) :: divs

      character (len=StrKIND) :: remainder_str
      character (len=StrKIND) :: temp_str

      type (mpas_time_type) :: ref_time
      type (mpas_timeinterval_type) :: int1, int2, remainder
      type (mpas_timeinterval_type) :: expected_remainder

      integer :: err_tmp

      ierr = 0

      call mpas_log_write('')
      call mpas_log_write('  Testing time intervals:')
      call mpas_log_write('      Reference time: ' // trim(ref_str))
      call mpas_log_write('      Interval 1: ' // trim(int1_str))
      call mpas_log_write('      Interval 2: ' // trim(int2_str))

      call mpas_set_time(ref_time, dateTimeString=ref_str, ierr=err_tmp)
      call mpas_set_timeinterval(int1, timeString=int1_str, ierr=err_tmp)
      call mpas_set_timeinterval(int2, timeString=int2_str, ierr=err_tmp)
      call mpas_set_timeinterval(expected_remainder, timeString=expected_remainder_str, ierr=err_tmp)

      call mpas_log_write('      -- Calling interval division')

      call mpas_interval_division(ref_time, int1, int2, divs, remainder)

      call mpas_get_timeinterval(remainder, startTimeIn=ref_time, timeString=remainder_str)

      call mpas_log_write('      Interval Division summary')
      write(temp_str,*) '         Divisions: ', divs
      call mpas_log_write(trim(temp_str))
      call mpas_log_write('          Remainder: ' // trim(remainder_str))
      call mpas_log_write('')

      if ( divs == expected_divs ) then
         call mpas_log_write('          Div Test: PASSED')
      else
         call mpas_log_write(' **       Div Test: FAILED', MPAS_LOG_ERR)
         ierr = 1
      end if

      if ( remainder == expected_remainder ) then
         call mpas_log_write('          Remainder Test: PASSED')
      else
         call mpas_log_write(' **       Remainder Test: FAILED', MPAS_LOG_ERR)
         ierr = 1
      end if


   end subroutine test_core_interval_test!}}}

end module test_core_timekeeping_tests
