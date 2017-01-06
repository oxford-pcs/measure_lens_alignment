# measure\_lens\_alignment

## The experiment 

It is probably useful to first talk a little bit about what is actually trying to be achieved with this software. The basic premise of the experiment is to ascertain the optical axes (OA, defined as a line connecting the two centres of curvature) and mechanical axes (MA, defined from the outer diameter of the lens itself) of the SWIFT lenses when mounted in their rings. Additionally, there is speculation that these axes are not consistent between rotations, suggesting that the lenses are moving in their mounts as the gravity vector changes. 

We therefore plan to measure the OA and MA for several orientations of the lens rings. We do this by first setting up a part coordinate system (PCS) to which all readings will be transformed into. This part coordinate system is defined by the front surface of the lens ring, the outer diameter of the lens ring and an M6 thru' hole on the left side of the mount (looking face on). A line between the outer diameter of the ring and the thru' hole defines the x/y axes, with the front surface defining the Z.

With this definition, one would expect that **FOR A SECURELY FIXED LENS** the x/y position of both the OA and MA (evaluated at Z=0) would draw a circle on an x/y position plot as the lens ring was rotated. Deviations from this circle indicate that the lens may be sagging. Also, hysteresis may cause the lens not to return to its original starting position on a full rotation.

## Setting up

This program processes **measurements** in a Quartis CMM database. Each **measurement** contains a group of **elements**; these are circles, spheres, lines etc. 

When taking a measurement, the Quartis CMM software has the option of attaching a PCS. The Quartis CMM database has two pertinent tables to record the information taken by the CMM: **\_tbCoordSys** and **\_tbElement**. The former table contains the matrix elements of a transform that can be used to convert the raw CMM table positions (e.g. **\_tbElement.elAct\***) into the desired PCS. The two tables can be joined by **\_tbElement.elMsRecNr --> \_tbCoordSys.csMsRecNr**, giving each element a PCS (many elements will have the same PCS). 

The configuration JSON file, config.json, provides the necessary contextual information for the program to proceed. Aside from general information regarding what sets of data to process, where the database resides and which tables to use (**PARAMS.lenses**, **PARAMS.db\_path**, **PARAMS.T\_COORDINATE\_SYSTEMS** and **PARAMS.T\_ELEMENTS** respectively), this file holds information regarding the different coordinate systems and what lens ring angles they were taken at. These are kept in the **COORDINATE\_SYSTEMS** array and are listed with an arbitrary **id** field, the angle of the lens ring at which the coordinate system was defined (**angle**) and the measurement record number it corresponds to in the database (**csMsRecNr**). 

**FOR A COORDINATE SYSTEM TO BE CONSIDERED, IT MUST BE ENTERED IN THIS SECTION OF THE CONFIGURATION FILE**.

To create a new dataset (a single dataset contains all rotations of the ring), a new entry in **DATA** is created, and a unique **id** given (which must be specified in **PARAMS.lenses** for it to be processed). A range of measurement ids is specified in **elMsRecNr\_range**, along with an array of elements (**elements**) the program can expect to find. Each entry in **elements** needs an arbitrary id, a db named element type (**elid** as in **\_tbElement.elId**, e.g. CIR\_1, PLN\_1, LIN\_1 etc.), a description (**desc**), an explicit element type (**type**, e.g CIRCLE, LINE, SPHERE but currently unused) and also a key (**key**). The key can be one of "MECH\_F", "MECH\_R", "OPT\_F", "OPT\_R" or "null" and defines if this particular measurement should be used in the calculation of the OA or MA (or not at all).

## How it works

Information regarding this is given as inline code in go.py.