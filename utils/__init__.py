"""Utility functions/exceptions sub-package"""
## \namespace geo.utils Pragmatic Utilities

import collections


## Thankyou, stackexchange.
class keydefaultdict(collections.defaultdict):
    """A default dict with function-of-key values."""

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        else:
            #pylint: disable=E1102
            ret = self[key] = self.default_factory(key)
            return ret
