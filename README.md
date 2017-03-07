# measure\_lens\_alignment

## Setting up the configuration file

The configuration JSON file "config.json" contains, amongst other things, the necessary contextual information defining a set of measurements. There are four key parameters that are contained with the 
configuration file: 

- PARAMS
- CONFIGURATIONS
- COORDINATE SYSTEMS, and
- DATA

**PARAMS** contains general information required for the program to run, e.g. where the database resides.

Each entry in **CONFIGURATIONS** defines a complete measurement set by specifying which datasets from **DATA** are to be used. Each has a unique identifier (**id**) and type (**type**).

For type "rotational\_analysis", additional fields are required for:

- an identifier of the dataset containing rotational data (**rotation\_data**)
- an identifier of the dataset containing error data (**error\_data**) 
- the element to be used for the front surface of the lens (**lens\_front\_elId**) 
- the element to be used for the rear surface of the lens (**lens\_rear\_elId**) 
- the element to be used for the mechanical front of the lens (**mount\_front\_elId**), 
- the element to be used for the mechanical rear of the lens (**mount\_rear\_elId**)
- the mount lens thickness (**mount\_ring\_thickness**)
- the PCS id to be used for both rotational and error data (**rotation_data_csId** and **error_data_csId**)
- the two angles to be used for hysteresis analysis (**hys\_idx\_1** and **hys\_idx\_2**), with indexes defined from the array of angles

**COORDINATE SYSTEM** entries add extra information to a database coordinate system entry. For example, some may have an "angle" field which specifies the angle of the lens mount at which the measurement 
was done and a csMsRecNr to allow it to be mapped to entries in the database. **FOR A COORDINATE SYSTEM TO BE CONSIDERED, IT MUST BE ENTERED IN THIS SECTION OF THE CONFIGURATION FILE**.

Measurement datasets are defined in **DATA**. Each dataset must have fields containing a unique identifier (**id**), a measurement record range (**elMsRecNr_range**) and a list of elements within each 
measurement. Each of these elements should have a unique identifier (**id**) field, an element ID corresponding to that defined in the database (**elId**) and a description field (**desc**).

## rotation_analysis.py

The basic premise of this routine is to calculate the optical axes (OA, defined as a line connecting the two centres of lens curvature) and mechanical axes (MA, defined from the outer diameter of the lens itself) 
of lenses when mounted in their rings. To check for stability, these calculations are done for different rotations (and thus gravity vectors) of the rings.

Obviously, the geometry of the measurements taken to ascertain the parameters that allow these two axes to be calculated will vary on the geometry of the setup. For SWIFT, we create a part coordinate system 
defined by the front surface of the lens ring, the outer diameter of the lens ring and an M6 thru' hole on the left side of the mount (looking face on). A line between the outer diameter 
of the ring and the thru' hole defines the x/y axes, with the front surface defining the z. Alongside this PCS, the program requires 4 elements to be measured: 

1. An element of type SPHERE representing the front surface of the lens
2. An element of type SPHERE representing the rear surface of the lens
3. An element of type CIRCLE representing the mechanical front diameter of the lens
4. An element of type CIRCLE representing the mechanical rear diameter of the lens

With this, the program can generate plots (-p2d) and information (-pi) regarding the calcuation of the mechanical axis (-ma) and optical axis (-oa) respectively by invoking e.g.

`python rotation_analysis.py -p2d -pi -ma -oa -l lens1 `


