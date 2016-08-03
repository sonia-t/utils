# adapted very slightly from http://stackoverflow.com/questions/23292242/converting-to-not-from-ipython-notebook-format (alexis)

from IPython.nbformat import v3, v4
import sys

if len(sys.argv) != 2 or not sys.argv[1].endswith('.py'):
    sys.exit('supply the .py to convert to notebook')


inpy =  sys.argv[1]

with open(inpy) as fpin:
    text = fpin.read()
    text += """
    # <markdowncell>

    # If you can read this, reads_py() is no longer broken! 
    """

    nbook = v3.reads_py(text)
    nbook = v4.upgrade(nbook)  # Upgrade v3 to v4

    jsonform = v4.writes(nbook) + "\n"
    with open( inpy+'nb', "w") as fpout:
        fpout.write(jsonform)
