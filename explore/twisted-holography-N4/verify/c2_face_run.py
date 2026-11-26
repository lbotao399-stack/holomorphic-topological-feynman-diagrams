import sys
from c2_counts import marked_face
for m in map(int, sys.argv[1:]):
    print(f'marked face formula m={m} (n={m-1}): {marked_face(m)}', flush=True)
