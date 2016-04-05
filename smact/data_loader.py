################################################################################
#  Copyright J. M. Sketon, D. W. Davies (2016)                                 #
#                                                                              #
#  This file is part of SMACT: smact.__init__ is free software: you can        #
#  redistribute it and/or modify it under the terms of the GNU General Public  #
#  License as published by the Free Software Foundation, either version 3 of   #
#  the License, or (at your option) any later version.                         #
#  This program is distributed in the hope that it will be useful, but WITHOUT #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   #
#  more details.                                                               #
#  You should have received a copy of the GNU General Public License along with#
#  this program.  If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

"""
This module handles the loading of external data used to initialise the core smact.Element and smact.Species classes.
It implements a transparent data-caching system to avoid a large amount of I/O when naively constructing several of these objects.
It also implements a switchable system to print verbose warning messages about possible missing data (mainly for debugging purposes).
"""


import csv;
import os;

from smact import data_directory;


# Module-level switch for printing "verbose" warning messages about missing data.

_PrintWarnings = False;

def ToggleWarnings(enable):
    """
    Toggles verbose warning messages on and off.
    ** In order to see any of the warnings, this function needs to be called _before_ the first call to the smact.Element() constructor. **
    
    Args:
        enable : print verbose warning messages.
    """
    
    global _PrintWarnings;
    
    _PrintWarnings = enable;


# Loader and cache for the element oxidation-state data.

_ElementOxidationStates = None;

def GetElementOxidationStates(symbol):
    """
    Retrieve a list of known oxidation states for an element.
    
    Args:
        symbol : the atomic symbol of the element to look up.
    
    Returns:
        A list of elements, or None if oxidation states for the element were not found in the external data.
    """
    
    global _ElementOxidationStates;
    
    if _ElementOxidationStates == None:
        _ElementOxidationStates = { };

        with open(os.path.join(data_directory, "oxidation_states.txt"), 'r') as file:
            for line in file:
                line = line.strip()

                if line[0] != '#':
                    items = line.split()

                    _ElementOxidationStates[items[0]] = [int(oxidationState) for oxidationState in items[1:]]

    if symbol in _ElementOxidationStates:
        return _ElementOxidationStates[symbol]
    else:
        if _PrintWarnings:
            print("WARNING: Oxidation states for element {0} not found.".format(symbol));
        
        return None;

# Loader and cache for the element crustal abundances.

_ElementCrustalAbundances = None;

def GetElementCrustalAbundance(symbol):
    """
    Retrieve the crustal abundance for an element.
    
    Args:
        symbol : the atomic symbol of the element to look up.
    
    Returns:
        The crustal abundance, or None if a value for the element was not found in the external data.
    """
    
    global _ElementCrustalAbundances;

    if _ElementCrustalAbundances == None:
        _ElementCrustalAbundances = { };

        with open(os.path.join(data_directory, "crustal_abundance.txt"), 'r') as file:
            for line in file:
                line = line.strip();

                if line[0] != '#':
                    items = line.split();

                    _ElementCrustalAbundances[items[0]] = float(items[1]);

    if symbol in _ElementCrustalAbundances:
        return _ElementCrustalAbundances[symbol];
    else:
        if _PrintWarnings:
            print("WARNING: Crustal-abundance data for element {0} not found.".format(symbol));
        
        return None;

# Loader and cache for the element HHI scores.

_ElementHHIs = None;

def GetElementHHIs(symbol):
    """
    Retrieve the HHI_R and HHI_p scores for an element.
    
    Args:
        symbol : the atomic symbol of the element to look up.
    
    Returns:
        A (HHI_p, HHI_R) tuple, or None if values for the elements were not found in the external data.
    """
    
    global _ElementHHIs;
    
    if _ElementHHIs == None:
        _ElementHHIs = { };

        with open(os.path.join(data_directory, "HHIs.txt"), 'r') as file:
            for line in file:
                items = line.strip();
                
                if line[0] != '#':
                    items = line.split();
                    

    if symbol in _ElementHHIs:
        return _ElementHHIs[symbol];
    else:
        if _PrintWarnings:
            print("WARNING: HHI data for element {0} not found.".format(symbol));
        
        return None;

# Loader and cache for the Open Babel-derived element data.

_ElementOpenBabelDerivedData = None;

def GetElementOpenBabelDerivedData(symbol):
    """
    Retrieve the Open Banel-derived data for an element.
    
    Args:
        symbol : the atomic symbol of the element to look up.
    
    Returns:
        A dictionary containing the data read from the Open Babel table, keyed with the column headings.
    """
    
    global _ElementOpenBabelDerivedData;
    
    if _ElementOpenBabelDerivedData == None:
        _ElementOpenBabelDerivedData = { };

        with open(os.path.join(data_directory, "element.txt"), 'r') as file:
            for line in file:
                line = line.strip();

                if line[0] != '#':
                    items = line.split();

                    key = items[1];

                    dataset = { };
                    
                    dataset['Number'] = int(items[0]);

                    areNeg = float(items[2]);
                    
                    # These values are not (presently) used by SMACT -> no need to print a warning.
                    
                    #if _PrintWarnings and areNeg == 0.0:
                    #    print("WARNING: Aldred-Rochow electronegativity for element {0} may be set to the default value of zero in the Open Babel-derived data.".format(symbol));

                    dataset['ARENeg'] = areNeg;

                    rCov = float(items[3]);

                    if _PrintWarnings and rCov == 1.6:
                        print("WARNING: Covalent radius for element {0} may be set to the default value of 1.6 in the Open Babel-derived data.".format(key));

                    dataset['RCov'] = float(items[3]);

                    dataset['RBO'] = float(items[4]);

                    rVDW = float(items[5]);

                    # These values are not (presently) used by SMACT -> no need to print a warning.

                    #if _PrintWarnings and rVDW == 2.0:
                    #    print("WARNING: Van der Waals raius for element {0} may be set to the default value of 2.0 in the Open Babel-derived data.".format(symbol));

                    dataset['RVdW'] = float(rVDW);
                    
                    maxBnd = int(items[6]);
                    
                    # These values are not (presently) used by SMACT -> no need to print a warning.

                    #if _PrintWarnings and maxBnd == 6:
                    #    print("WARNING: Maximum bond valence for element {0} may be set to the default value of 6 in the Open Babel-derived data.".format(symbol));

                    dataset['MaxBnd'] = maxBnd;

                    dataset['Mass'] = float(items[7]);
                    
                    elNeg = float(items[8]);

                    if _PrintWarnings and elNeg == 0.0:
                        print("WARNING: Pauling electronegativity for {0} may be set to the default of zero in the Open Babel-derived data.".format(key));
                    
                    dataset['ElNeg.'] = elNeg;

                    ionization = float(items[9]);

                    if _PrintWarnings and ionization == 0.0:
                        print("WARNING: Ionisation potential for {0} may be set to the default of zero in the Open Babel-derived data.".format(key));

                    dataset['Ionization'] = ionization;
                    
                    elAffinity = float(items[10]);

                    if _PrintWarnings and elAffinity == 0.0:
                        print("WARNING: Electron affinity for {0} may be set to the default of zero in the Open Babel-derived data.".format(key));
                    
                    dataset['ElAffinity'] = elAffinity;

                    dataset['RGB'] = (float(items[11]), float(items[12]), float(items[13]));
                    dataset['Name'] = items[14];

                    _ElementOpenBabelDerivedData[items[1]] = dataset;

    if symbol in _ElementOpenBabelDerivedData:
        return _ElementOpenBabelDerivedData[symbol];
    else:
        if _PrintWarnings:
            print("WARNING: Open Babel-derived element data for {0} not found.".format(symbol));
        
        return None;

# Loader and cache for the element eigenvalues.

_ElementEigenvalues = None;

def GetElementEigenvalue(symbol):
    """
    Retrieve the eigenvalue for an element.
    
    Args:
        symbol : the atomic symbol of the element to look up.
    
    Returns:
        The eigenvalue, or None if an eigenvalue was not found in the external data.
    """
    
    global _ElementEigenvalues;
    
    if _ElementEigenvalues == None:
        _ElementEigenvalues = { };
        
        with open(os.path.join(data_directory, "Eigenvalues.csv"),'r') as file:
            reader = csv.reader(file);

            # Skip the first row (headers).
            
            next(reader);

            for row in reader:
                _ElementEigenvalues[row[0]] = float(row[1]);
    
    if symbol in _ElementEigenvalues:
        return _ElementEigenvalues[symbol];
    else:
        if _PrintWarnings:
            print("WARNING: Eigenvalue data for element {0} not found.".format(symbol));
        
        return None;

# Loader and cache for the element s eigenvalues.

_ElementSEigenvalues = None;

def GetElementSEigenvalue(symbol):
    """
    Retrieve the s eigenvalue for an element.
    
    Args:
        symbol : the atomic symbol of the element to look up.
    
    Returns:
        The s eigenvalue, or None if an s eigenvalue was not found in the external data.
    """
    
    global _ElementSEigenvalues;
    
    if _ElementSEigenvalues == None:
        _ElementSEigenvalues = { };
        
        with open(os.path.join(data_directory, "Eigenvalues_s.csv"), 'rU') as file:
            reader = csv.reader(file);

            # Skip the first row (headers).
            
            next(reader);

            for row in reader:
                _ElementEigenvalues[row[0]] = float(row[1]);

    if symbol in _ElementSEigenvalues:
        return _ElementSEigenvalues[symbol];
    else:
        if _PrintWarnings:
            print("WARNING: s-eigenvalue data for element {0} not found.".format(symbol));
        
        return None;

# Loader and cache for the element Shannon radii datasets.

_ElementShannonRadiiData = None;

def GetElementShannonRadiusData(symbol):
    """
    Retrieve Shannon radii for known oxidation states/coordination environments of an element.
    
    Args:
        symbol : the atomic symbol of the element to look up.
    
    Returns:
        A list of Shannon radii datasets, or None if the element was not found among the external data.
    """
    
    global _ElementShannonRadiiData;
    
    if _ElementShannonRadiiData == None:
        _ElementShannonRadiiData = { };
        
        with open(os.path.join(data_directory, "shannon_radii.csv"), 'rU') as file:
            reader = csv.reader(file);
            
            # Skip the first row (headers).
            
            next(reader);
        
            for row in reader:
                # For the shannon radii, there are multiple datasets for different element/oxidation-state/coordination combinations.
                
                key = row[0];
                
                dataset = {
                    'charge' : int(row[1]),
                    'coordination' : row[2],
                    'crystal_radius' : float(row[3]),
                    'ionic_radius' : float(row[4]),
                    'comment' : row[5]
                    };
                
                if key in _ElementShannonRadiiData:
                    _ElementShannonRadiiData[key].append(dataset);
                else:
                    _ElementShannonRadiiData[key] = [dataset];
    
    if symbol in _ElementShannonRadiiData:
        return _ElementShannonRadiiData[symbol];
    else:
        if _PrintWarnings:
            print("WARNING: Shannon-radius data for element {0} not found.".format(symbol));
        
        return None;

# Loader and cache for the element solid-state energy (SSE) datasets.

_ElementSSEData = None;

def GetElementSSEData(symbol):
    """
    Retrieve the solid-state energy (SSE) data for an element.
    
    Args:
        symbol : the atomic symbol of the element to look up.
    
    Returns:
        A dictionary containing the SSE dataset for the element, or None if the element was not found among the external data.
    """
    
    global _ElementSSEData;
    
    if _ElementSSEData == None:
        _ElementSSEData = { };
        
        with open(os.path.join(data_directory, "SSE.csv"), 'rU') as file:
            reader = csv.reader(file)
            
            for row in reader:
                dataset = {
                    'AtomicNumber' : int(row[1]),
                    'SolidStateEnergy' : float(row[2]),
                    
                    #TODO: Someone must know what these fields are -- and they might be useful!
                    
                    'Unknown1' : float(row[3]),
                    'Unknown2' : float(row[4]),
                    'Unknown3' : float(row[5]),
                    'Unknown4' : float(row[6])
                    };
                
                _ElementSSEData[row[0]] = dataset;

    if symbol in _ElementSSEData:
        return _ElementSSEData[symbol];
    else:
        if _PrintWarnings:
            print("WARNING: Solid-state energy data for element {0} not found.".format(symbol));
        
        return None;
