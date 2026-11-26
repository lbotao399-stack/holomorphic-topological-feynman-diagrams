import sys
p = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code'
if p not in sys.path:
    sys.path.insert(0, p)
