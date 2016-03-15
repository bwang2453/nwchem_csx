from nwchem import *

rtdb_open("h2o.db","old")

rtdb_put("marat",1)
rtdb_put("geometry:geometry:ncenter1",3)
n1=rtdb_get("geometry:geometry:ncenter1")
#print n1
#rtdb_print(1)
ncenter=rtdb_get("geometry:geometry:ncenter")
print ncenter
coords = rtdb_get('geometry:' + 'geometry' + ':coords')
print coords
#
atomNum = rtdb_get('geometry:' + 'geometry' + ':ncenter')
atmName = rtdb_get('geometry:' + 'geometry' + ':tags')
charge = rtdb_get('geometry:' + 'geometry' + ':charges')
#
print "Number of Atoms =",atomNum
print "Atom Names =",atmName
print "Atom Charges = ", charge
