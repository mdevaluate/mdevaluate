from libc.math cimport floor, ceil, fabs, log10

def is_log_step(long long step, long long  per_decade):

    # The first step will never be written
    if step == 0:
        return True

    if step == 1:
        return True


    cdef double log_step = log10(step)
    cdef double log_next = log10(step+1)
    cdef double log_prev = log10(step-1)
    cdef long long decade = int(log10(step))

    # Idea: Is there a number x = decade + n / per_decade
    # with 0 <= n < per_decade where
    # | log_step - x | <= | log_prev - x |
    # | log_step - x | <  | log_next - x |

    # Calculate the closest numbers

    cdef double n_smaller = floor((log_step - decade) * per_decade);
    cdef double n_larger = ceil((log_step - decade) * per_decade);

    cdef double prev_dist = fabs(log_prev-decade-n_smaller/per_decade);
    cdef double dist_to_smaller = fabs(log_step-decade-n_smaller/per_decade);

    cdef double next_dist = fabs(log_next-decade-n_larger/per_decade);
    cdef double dist_to_larger = fabs(log_step-decade-n_larger/per_decade);

    if dist_to_larger < next_dist or dist_to_smaller <= prev_dist:
        return True
    else:
        return False
