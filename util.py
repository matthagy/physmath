

def memorize(func):
    cache = {}
    def wrap(*args, **kwds):
        key = (args,() if not kwds else tuple(sorted(kwds.iteritems())))
        try:
            return cache[key]
        except KeyError:
            v = cache[key] = func(*args, **kwds)
            return v
    wrap.clear = cache.clear
    return wrap
