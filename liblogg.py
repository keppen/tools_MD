# import logging
import time


def timing(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time

        hours = int(elapsed_time / 3600)
        minutes = int((elapsed_time % 3600) / 60)
        seconds = int(elapsed_time % 60)
        miliseconds = int((elapsed_time - int(elapsed_time)) * 1e3)

        print(f"Function '{func.__name__}' took {hours:02d}:{minutes:02d}:{seconds:02d}.{miliseconds:3} seconds to execute.")
        return result
    return wrapper
