'''Utilitiy to reload project
'''

from __future__ import absolute_import

def fixup(force=False):
    import sys
    import os
    import time
    from types import ModuleType

    load_order = '''
    sigfig
    algebra
    dne
    ratio
    units
    physnum
    types
    layout
    convert
    calculate
    '''.split()

    def main():
        now = time.time()
        global last_load_time
        try:
            last_load_time
        except NameError:
            last_load_time = None
        reloaded = maybe_nuke_jamenson(force)
        for name in load_order:
            fullname = 'physmath.'+name
            if fullname not in sys.modules:
                continue
            if not reloaded:
                tm = get_module_mtime(fullname)
                #print fullname, tm - last_load_time
                if tm < last_load_time:
                    continue
            msg('reloading %s', name)
            reload(sys.modules[fullname])
            reloaded = True
        last_load_time = now

        #import browsemath.client
        #reload(browsemath.client)

    def msg(msg, *args):
        print >>sys.stderr, msg%args if args else msg

    def maybe_nuke_jamenson(force):
        '''Dependencies here are tricky, find if anything has changed and nuke everything
           if needed
        '''
        global last_load_time
        names = list(name for name in sys.modules if name.startswith('jamenson'))
        if not force and last_load_time is not None:
            mtime = max(filter(None, map(get_module_mtime, names)))
            if mtime < last_load_time:
                return False
        msg('nuking jamenson')
        map(sys.modules.pop, names)
        return True

    def get_module_mtime(name):
        mod = sys.modules[name]
        if not isinstance(mod, ModuleType):
            modname = getattr(type(mod), '__name__', None)
            if modname and modname != name and modname in sys.modules:
                return get_module_mtime(modname)
            return 0
        path = mod.__file__
        if path.endswith('pyc'):
            path = path[:-1:]
        while os.path.islink(path):
            path = os.readlink(path)
        #print path
        try:
            st = os.stat(path)
        except OSError:
            return 0xffffffff
        return st.st_mtime
    main()

fixup()
