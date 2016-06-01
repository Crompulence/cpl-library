import os
import shutil as sh

dirlist = os.listdir('./')

deletelist = []
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

answer = raw_input('Proceed? [y]/n: ')
if (answer == 'y' or answer == 'Y'):
    for f in deletelist:
        sh.rmtree(f)
else:
    quit('Cancelled deletion.')
