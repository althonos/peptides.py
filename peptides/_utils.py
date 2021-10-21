"""Internal utility functions for the `peptides` package.
"""

def descriptor(prefix):

    def decorator(func):
        func.descriptor = True
        func.prefix = prefix
        return func

    return decorator
