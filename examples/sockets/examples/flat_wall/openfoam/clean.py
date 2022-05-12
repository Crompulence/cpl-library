import os
import shutil as sh

dirlist = os.listdir('./')

MESH_DIR = './constant/polyMesh'
MESH_DEL_FILES = ['boundary', 'faces', 'neighbour', 'owner', 'points']
deletelist = [os.path.join(MESH_DIR, file) for file in MESH_DEL_FILES]
keeplist = []

for f in dirlist:

    if (f == '0'):
        continue

    # Delete if filename can be turned into a float and is not 0
    try:

        testfloat = float(f)
        deletelist.append(f)

    except ValueError:

        if ('processor' in f):
            deletelist.append(f)

keeplist = [i for i in dirlist if i not in deletelist]


print('Keeping the following files:')
print(keeplist)

print('Deleting the following files: ')
print(deletelist)

answer = input('Proceed? [y]/n: ')
if (answer == 'y' or answer == 'Y'):
    for f in deletelist:
        try:
            sh.rmtree(f)
        except:
            try:
                os.remove(f)
            except:
                pass
else:
    quit('Cancelled deletion.')
