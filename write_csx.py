import fortranfile
import numpy as np
import csx1_api as api
import sys
import math
import struct
from nwchem import *

dbFile = sys.argv[1]
rtdb_open(dbFile, 'old')
#determine the calculation type
hasFreq = False
hasProp = False
calcType = rtdb_get('task:' + 'theory')
try:
    calcOper = rtdb_get('task:' + 'operation')
    if calcOper == 'freq':
        hasFreq = True
    if calcOper == 'property':
        hasProp = True
except NwchemError:
    pass
fileRoot = rtdb_get('file_prefix')
csxfile = open( fileRoot + '.csx', 'w')

#get the coordinates
coords = rtdb_get('geometry:' + 'geometry' + ':coords')
autoang = 1.889725989
atomNum = rtdb_get('geometry:' + 'geometry' + ':ncenter')
atmName = rtdb_get('geometry:' + 'geometry' + ':tags')
charge = rtdb_get('geometry:' + 'geometry' + ':charges')
atmMass = rtdb_get('geometry:' + 'geometry' + ':masses')
try:
    xcType = rtdb_get('dft:' + 'xc_spec')[0]
except NwchemError:
    xcType = 'slater vwn_5'
try: 
    molCharge = rtdb_get('charge')
except NwchemError:
    molCharge = 0.0
molCharge = int(molCharge)
if calcType == 'dft' or calcType == 'tddft':
    molMulti = rtdb_get('dft:' + 'nopen')+1
else: 
    molMulti = rtdb_get('scf:' + 'nopen')+1

wfnRestricted = False
if molMulti == 1:
    molSpin = 'RHF'
    wfnRestricted = True
else:
    molSpin = 'UHF'
try: 
    # If this fails then the basis set was probably specified with a star (*)
    basName = rtdb_get('basis:' + 'ao basis' + ':bs_stdname')[0].upper()
except NwchemError:
    # If this also fails then no basis set was specified???
    basName = rtdb_get('basis:' + 'ao basis' + ':star bas type').upper()
#Wavefunction
filename = fileRoot + '.movecs'
mofile = fortranfile.FortranFile(filename)
x = mofile.readString()
y = mofile.readString().strip()
z = mofile.readInts()
a = mofile.readString()
b = mofile.readInts()
c = mofile.readString()
d = mofile.readInts()
if wfnRestricted:
    orbNum = mofile.readInts()[0]
    nbf = mofile.readInts()[0]
    orbOcc = mofile.readReals()
    orbOccString = ' '.join(str(x) for x in sorted(orbOcc,reverse=True))
    orbEne = mofile.readReals('d')
    orbEString = ' '.join(str(x) for x in orbEne)
    orbCaString = []
    for iorb in range(orbNum):
        orbCa = mofile.readReals('d')
        orbCaString.append(' '.join(str(x) for x in orbCa))
    wfn1 = api.waveFunctionType(orbitalCount=orbNum,orbitalOccupancies=orbOccString)
    orbe1 = api.stringArrayType(unit='cs:hartree')
    orbe1.set_valueOf_(orbEString)
    orbs1 = api.orbitalsType()
    for iorb in range(orbNum):
            orbt = orbCaString[iorb]
            orb1 = api.stringArrayType(id=iorb+1)
            orb1.set_valueOf_(orbt)
            orbs1.add_orbital(orb1)
    wfn1.set_orbitals(orbs1)
    wfn1.set_orbitalEnergies(orbe1)
else:
    orbNum = mofile.readInts()[0]
    nbf = mofile.readInts()[0]
    orbOccCa = mofile.readReals('d')
    orbOccCaString = ' '.join(str(x) for x in sorted(orbOccCa,reverse=True))
    orbCaEne = mofile.readReals('d')
    orbCaEString = ' '.join(str(x) for x in orbCaEne)
    orbCaString = []
    for iorb in range(orbNum):
        orbCa = mofile.readReals('d')
        orbCaString.append(' '.join(str(x) for x in orbCa))
    orbOccCb = mofile.readReals('d')
    orbOccCbString = ' '.join(str(x) for x in sorted(orbOccCb,reverse=True))
    orbCbEne = mofile.readReals('d')
    orbCbEString = ' '.join(str(x) for x in orbCbEne)
    orbCbString = []
    for iorb in range(orbNum):
        orbCb = mofile.readReals('d')
        orbCbString.append(' '.join(str(x) for x in orbCb))
    wfn1 = api.waveFunctionType(orbitalCount=orbNum)
    orbe1 = api.stringArrayType(unit='cs:hartree')
    orbe1.set_valueOf_(orbCaEString)
    orbs1 = api.orbitalsType()
    for iorb in range(orbNum):
            orbt = orbCaString[iorb]
            orb1 = api.stringArrayType(id=iorb+1)
            orb1.set_valueOf_(orbt)
            orbs1.add_orbital(orb1)
    wfn1.set_alphaOrbitals(orbs1)
    wfn1.set_alphaOrbitalEnergies(orbe1)
    wfn1.set_alphaOrbitalOccupancies(orbOccCaString)
    orbe2 = api.stringArrayType(unit='cs:hartree')
    orbe2.set_valueOf_(orbCbEString)
    orbs2 = api.orbitalsType()
    for iorb in range(orbNum):
            orbt = orbCbString[iorb]
            orb2 = api.stringArrayType(id=iorb+1)
            orb2.set_valueOf_(orbt)
            orbs2.add_orbital(orb2)
    wfn1.set_betaOrbitals(orbs2)
    wfn1.set_betaOrbitalEnergies(orbe2)
    wfn1.set_betaOrbitalOccupancies(orbOccCbString)

#dipole moment information
if hasProp:
    dipole = rtdb_get('dft:' + 'dipole')
    au2db = 2.541766
    molDipoleX = dipole[0]*au2db
    molDipoleY = dipole[1]*au2db
    molDipoleZ = dipole[2]*au2db
    molDipoleTot = math.sqrt(molDipoleX*molDipoleX+molDipoleY*molDipoleY+molDipoleZ*molDipoleZ)
    prop1 = api.propertiesType()
    sprop1 = api.propertyType(name='dipoleMomentX',unit='cs:debye')
    sprop1.set_valueOf_(molDipoleX)
    sprop2 = api.propertyType(name='dipoleMomentY',unit='cs:debye')
    sprop2.set_valueOf_(molDipoleY)
    sprop3 = api.propertyType(name='dipoleMomentZ',unit='cs:debye')
    sprop3.set_valueOf_(molDipoleZ)
    sprop4 = api.propertyType(name='dipoleMomentAverage',unit='cs:debye')
    sprop4.set_valueOf_(molDipoleTot)
    prop1.add_systemProperty(sprop1)
    prop1.add_systemProperty(sprop2)
    prop1.add_systemProperty(sprop3)
    prop1.add_systemProperty(sprop4)

#vibration information
if hasFreq:
    normFileName = fileRoot + '.nmode'
    normFile = open(normFileName, 'rb')
    normNum = int(struct.unpack("@d", normFile.read(8))[0])
    molFreqNum = int(struct.unpack("@d", normFile.read(8))[0])
#       vibFreqs = rtdb_get('vib:'+'frequencies')
    vibFreqs = []
    for ifrq in range(molFreqNum):
        vibFreqs.append(struct.unpack("@d", normFile.read(8))[0])
    frqString = ' '.join(str(x) for x in vibFreqs)
    vibInts = rtdb_get('vib:'+'intensities')
    intString = ' '.join(str(x) for x in vibInts)
    vib1 = api.vibAnalysisType(vibrationCount=molFreqNum)
    freq1 = api.stringArrayType(unit="cs:cm-1")
    freq1.set_valueOf_(frqString)
    vib1.set_frequencies(freq1)
    irint1 = api.stringArrayType()
    irint1.set_valueOf_(intString)
    vib1.set_irIntensities(irint1)
    norms1 = api.normalModesType()
    normMdString = []
    for ifrq in range(molFreqNum):
        normM = []
        for inorm in range(normNum):
            normM.append(struct.unpack("@d", normFile.read(8))[0])
        normMdString.append(' '.join(str(x) for x in normM))
        norm1 = api.normalModeType(id=ifrq+1)
        norm1.set_valueOf_(normMdString[ifrq])
        norms1.add_normalMode(norm1)
    normFile.close()
    vib1.set_normalModes(norms1)

#Start to generate CSX elements
cs1 = api.csType(version='1.0')

#molecular publication section
mp1 = api.mpType(title=fileRoot, \
        abstract='default abstract', \
        publisher='default publisher', \
        status='default status', \
        category=1, \
        visibility=0, \
        tags='NWCHEM', \
        key=1 )

source1 = api.sourcePackageType(name='NWCHEM', version='6.5')
mp1.set_sourcePackage(source1)
ath1 = api.authorType(creator='bwang', \
        type_='cs:corresponding', \
        organization='default organization', \
        email='bwang@chemicalsemantics.com')
mp1.add_author(ath1)
cs1.set_molecularPublication(mp1)

#molecular system section
ms1 = api.msType(systemCharge=molCharge, \
       systemMultiplicity=molMulti)
temp1 = api.dataWithUnitsType(unit='cs:kelvin')
temp1.set_valueOf_(0.0)
ms1.set_systemTemperature(temp1)
mol1 = api.moleculeType(id='m1',atomCount=atomNum)
#obmol1 = openbabel.OBMol()
#etb = openbabel.OBElementTable()
for iatm in range(atomNum):
    #   xCoord = float(data.atomcoords[iatm,0])
    xCoord = coords[3*iatm]/autoang
    yCoord = coords[3*iatm+1]/autoang
    zCoord = coords[3*iatm+2]/autoang
    xCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
    yCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
    zCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
    xCoord1.set_valueOf_(xCoord)
    yCoord1.set_valueOf_(yCoord)
    zCoord1.set_valueOf_(zCoord)
    atm = api.atomType(id='a'+str(iatm+1), \
            elementSymbol=atmName[iatm], \
            atomMass=atmMass[iatm], \
            xCoord3D=xCoord1, \
            yCoord3D=yCoord1, \
            zCoord3D=zCoord1, \
            basisSet='cs:' + basName, \
            calculatedAtomCharge=0, \
            formalAtomCharge=0)
    mol1.add_atom(atm)
ms1.add_molecule(mol1)
cs1.set_molecularSystem(ms1)

#molCalculation section
sd_wfn_method = ['scf', 'dft', 'mp2', 'mp3', 'mp4']
md_wfn_method = ['ccd', 'ccsd', 'ccsdt']
mc1 = api.mcType()
qm1 = api.qmCalcType()
srs1 = api.srsMethodType()
if calcType in sd_wfn_method:
    sdm1 = api.srssdMethodType()
    #DFT
    if (calcType == 'dft') or (calcType == 'tddft'):
        molEE = rtdb_get('dft:'+'energy')
        try:
            xcType = rtdb_get('dft:'+'xc_spec')
    #        xcType = rtdb_get('dft:'+'xc_spec')[0] + rtdb_get('dft:'+'xc_spec')[1]
        except NwchemError:
            xcType = 'slater vwn_5'
        theo1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                basisSet='bse:'+basName, dftFunctional='cs:'+xcType)
        ene1 = api.energiesType(unit='cs:hartree')
        ee_ene1 = api.energyType(type_='cs:totalPotential')
        ee_ene1.set_valueOf_(float(molEE))
        ene1.add_energy(ee_ene1)
        theo1.set_energies(ene1)
        theo1.set_waveFunction(wfn1)
        if hasProp:
            theo1.set_properties(prop1)
        if hasFreq:
            theo1.set_vibrationalAnalysis(vib1)
        sdm1.set_dft(theo1)
    #SCF
    elif (calcType == 'scf'):
        molEE = rtdb_get('scf:'+'energy')
        theo1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                basisSet='bse:'+basName)
        ene1 = api.energiesType(unit='cs:hartree')
        ee_ene1 = api.energyType(type_='cs:totalPotential')
        ee_ene1.set_valueOf_(float(molEE))
        ene1.add_energy(ee_ene1)
        theo1.set_energies(ene1)
        theo1.set_waveFunction(wfn1)
        if hasProp:
            theo1.set_properties(prop1)
        if hasFreq:
            theo1.set_vibrationalAnalysis(vib1)
        sdm1.set_abinitioScf(theo1)
    #MP2
    elif (calcType == 'mp2'):
        molEE = rtdb_get('mp2:'+'energy')
        molCE = rtdb_get('mp2:'+'correlation energy')
        theo1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                basisSet='bse:'+basName)
        ene1 = api.energiesType(unit='cs:hartree')
        ee_ene1 = api.energyType(type_='cs:totalPotential')
        ee_ene1.set_valueOf_(float(molEE))
        ce_ene1 = api.energyType(type_='cs:correlation')
        ce_ene1.set_valueOf_(float(molCE))
        ene1.add_energy(ee_ene1)
        ene1.add_energy(ce_ene1)
        theo1.set_energies(ene1)
        theo1.set_waveFunction(wfn1)
        if hasProp:
            theo1.set_properties(prop1)
        if hasFreq:
            theo1.set_vibrationalAnalysis(vib1)
        sdm1.set_mp2(theo1)
    else:
        print ('The current CSX does not support this method')
    srs1.set_singleDeterminant(sdm1)

if calcType in md_wfn_method:
    mdm1 = api.srsmdMethodType()
    if (calcType == 'ccsd'):
        molEE = rtdb_get('ccsd:'+'energy')
        molCE = rtdb_get('ccsd:'+'ccsd correlation energy')
        theo1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                basisSet='bse:'+basName)
        ene1 = api.energiesType(unit='cs:hartree')
        ee_ene1 = api.energyType(type_='cs:totalPotential')
        ee_ene1.set_valueOf_(float(molEE))
        ce_ene1 = api.energyType(type_='cs:correlation')
        ce_ene1.set_valueOf_(float(molCE))
        ene1.add_energy(ee_ene1)
        ene1.add_energy(ce_ene1)
        theo1.set_energies(ene1)
        theo1.set_waveFunction(wfn1)
        if hasProp:
            theo1.set_properties(prop1)
        if hasFreq:
            theo1.set_vibrationalAnalysis(vib1)
        mdm1.set_ccsd(theo1)
    elif (calcType == 'ccsdt'):
        molEE = rtdb_get('ccsdt:'+'energy')
        molCE = rtdb_get('ccsdt:'+'correlation energy')
        theo1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                basisSet='bse:'+basName)
        ene1 = api.energiesType(unit='cs:hartree')
        ee_ene1 = api.energyType(type_='cs:totalPotential')
        ee_ene1.set_valueOf_(float(molEE))
        ce_ene1 = api.energyType(type_='cs:correlation')
        ce_ene1.set_valueOf_(float(molCE))
        ene1.add_energy(ee_ene1)
        ene1.add_energy(ce_ene1)
        theo1.set_energies(ene1)
        theo1.set_waveFunction(wfn1)
        if hasProp:
            theo1.set_properties(prop1)
        if hasFreq:
            theo1.set_vibrationalAnalysis(vib1)
        mdm1.set_ccsd_t(theo1)
    else:
        print ('The current CSX does not support this method')
    srs1.set_multipleDeterminant(mdm1)

qm1.set_singleReferenceState(srs1)
mc1.set_quantumMechanics(qm1)
cs1.set_molecularCalculation(mc1)

csxfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
cs1.export(csxfile,0)
csxfile.close()
