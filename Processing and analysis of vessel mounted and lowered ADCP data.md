# Processing and analysis of vessel mounted and lowered ADCP datasets

This document will guide you through the steps necessary to process and  analyse both vessel mounted (VM-ADCP) and lowered acoustic doppler current profilers (L-ADCP). Both are complex data sources that contain a lot of noise and potential biases, but dont despair, with todays programms and code packages anyone can work with this data and produce meaningfull results.

We used the following data sources:

from the field:
- VM-ADCP: RDI Workhorse 75 kHz and 150 kHz (http://www.teledynemarine.com/workhorse-mariner-adcp?BrandID=16)
- L-ADCP: RDI Workhorse Sentinel 300kHz (http://www.teledynemarine.com/workhorse-sentinel-adcp?BrandID=16)
- Sea Bird 911plus CTD profiles (http://www.seabird.com/sbe911plus-ctd)

from the web:
- International Bathymetric Chart of the Arctic Ocean IBCAO (https://www.ngdc.noaa.gov/mgg/bathymetry/arctic/arctic.html)
- Global coastline dataset (https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/)

We used the following software:

Data collection:
- VMDAS ADCP data collection software (http://www.teledynemarine.com/rdi/support#)
- BB Talk terminal software (should also work with other terminal programmes, http://www.teledynemarine.com/rdi/support#)  
- Sea Bird CDT data collection suite (Seaterm, Seasave, SBE Data processing, http://www.seabird.com/sbe911plus-ctd)

Post processing:
- Oracle VM VirtualBox 5.0 running a UNIX machine with the CODAS ADCP processing suite (https://currents.soest.hawaii.edu/docs/adcp_doc/codas_doc/)

Analysis:
- Matlab 2016a, incl. the mapping and statistics toolbox and the following additional packages
- m_map (https://www.eoas.ubc.ca/~rich/map.html)
- seawater liberary (http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm)
- objective mapping (http://mooring.ucsd.edu/software/matlab/doc/toolbox/datafun/objmap.html)
- L-ADCP processing suite LDEO (http://www.ldeo.columbia.edu/~ant/LADCP)

Here is a visualization of VM-ADCP and L-ADCP current profiles:

![Here is a visualization of VM-ADCP and L-ADCP current profiles](ladcp_map_sections_3d_2014_with_vmadcp.png)

# Fieldwork

This section explains how ADCP data was gathered during each expedition.

## VM-ADCP

The VM-ADCP is usualy mounted in the ships hull or a retractable keel, and connected to an onboard PC. To record data one of two proramms is usually used: The Windows programm VMDAS (Vesse mounted data aquisition system) from RDI instruments or the open source suite UHDAS from the University of Hawaii.

Make sure that the PC is running with the exact UTC time. If you can, synchronize the PC clock with a GPS device, the NMEA navigation data stream or the internet. Having the correct timestamps in the ADCP data is important for the post-processing (combining the L-ADCP, VM-ADCP and CTD data).

<img src="vmdas.png">

Here we used VMDAS to collect VM-ADCP data. After starting up the PC and ensuring the PC is properly connected to the ADCP and the ships navigation channels, one can open the VMDAS and select the button: File--> Collect Data. This will open some empty diagrams. The next step is to edit the data options under: Options --> Edit data options. In the first section "Communications" the ADCP and NMEA (navigation) data channels need to be set to the right ports, but the ships technician likely did this already. In the next section "ADCP setup" you need to select the bin number and size. This depends on the frequency of the ADCP. We used 100 bins with a binsize of 8 m and 8 ms paceing between the bins. The ADCP also requires a specific set up file, which states details of the ADCPs mounting and ping-ing patters. This is usually avialable from the ships technician or stored on the data collection PC. Here is an example from RV Helmer Hanssen:

```
-----------------------------------------------------------------------------\
; ADCP Command File for use with VmDas software.
;
; ADCP type:     75 Khz Ocean Surveyor
; Setup name:    default
; Setup type:    Low resolution, long range profile(narrowband)
;
; NOTE:  Any line beginning with a semicolon in the first
;        column is treated as a comment and is ignored by
;        the VmDas software.
;
; NOTE:  This file is best viewed with a fixed-point font (e.g. courier).
; Modified Last: 12August2003
;----------------------------------------------------------------------------/

; Restore factory default settings in the ADCP
cr1

; set the data collection baud rate to 38400 bps,
; no parity, one stop bit, 8 data bits
; NOTE:  VmDas sends baud rate change command after all other commands in
; this file, so that it is not made permanent by a CK command.
cb611

;
;
;	CX Trigger input, default is CX0,0 no trigging
;      
;	CX1,0 waits for a positive trigging signal
CX1,0
;
;

; Set for narrowband single-ping profile mode (NP), hundred (NN) 8 meter bins (NS),
; 8 meter blanking distance (NF)
WP0
NN100
NP00001
NS0800
NF0800

; Disable single-ping bottom track (BP),
; Set maximum bottom search depth to 1000 meters (BX)
BP000
BX10000

; output velocity, correlation, echo intensity, percent good
ND111110000

; One second between bottom and water pings
TP000000

; Zero seconds between ensembles
; Since VmDas uses manual pinging, TE is ignored by the ADCP.
; You must set the time between ensemble in the VmDas Communication options
TE00000000

; Not set to calculate speed-of-sound, no depth sensor, external synchro heading
; sensor, no pitch or roll being used, no salinity sensor, use internal transducer
; temperature sensor
EZ0000001
;
; Sets the speed of sound to 1500m/s
;
EC1500
; Output beam data (rotations are done in software)
EX00000

; Set transducer misalignment (hundredths of degrees)
EA04530

; Set transducer depth (decimeters)
ED0042

; Set Salinity (ppt)
ES35

; save this setup to non-volatile memory in the ADCP
CK
```
The next section "Recording" includes setting the naming and storage options for the data. The section "NAV" lets you select the correct NMEA channels to co-record with the ADCP data. Storing the correct NMEA data is important for the post-processing, since we want to remove biases that stem from the ships movement. The other sections are less important and you can now press the OK button and start recording ADCP data by pressing the play button in the upper left corner.

During the expedition  make sure that the ADCP does not interfere with other acoustic sensors, and adjust the pining intervals if necessary.

At the end of the expedition, open the data collection folder and copy the entire folder onto you post processing PC. We need especially the ENR, N1R and N2R files.

## L-ADCP

The lowered ADCP should be firmly installed onto the CTD frame. Insert a battery with sufficient voltage (>30V) in the ADCP and connect the ADCP to a Windows PC via a serial port.

Make sure that the PC is running with the exact UTC time. If you can, synchronize the PC clock with a GPS device, the NMEA navigation data stream or the internet. Having the correct timestamps in the ADCP data is important for the post-processing (combining the L-ADCP, VM-ADCP and CTD data).

Following are instructions for launch and recovery.

Launch:

- Open the BB talk program, select the port the ADCP is connected to and open a terminal window
- wake up the ADCP by pressing <end> on the keyboard
- Erase old data on the instrument: type `RE ErAsE` and confirm by pressing <ENTER> (Do this if you are sure the data was saved before!)
- From menu: Transfer --> PC time to ADCP
- Start the ADCP by transferring the following ADCP settings file onto the ADCP: Press <F2> . Choose `ladcp_setting.txt` and wait until script ends
- Disconnect cables and launch CTD

create `ladcp_setting.txt` by saving the following code:
```
$P *************************************************************************
$P *******  LADCP Deployment with one ADCP.  Looking down **********
$P *************************************************************************
; Send ADCP a BREAK
$B
; Wait for command prompt (sent after each command)
$W62
; Display real time clock setting
tt?
$W62
; Set to factory defaults
CR1
$W62
; use WM15 for firmware 16.3
; activates LADCP mode (BT from WT pings)
WM15
; Flow control (Record data internally):
; - automatic ensemble cycling (next ens when ready)
; - automatic ping cycling (ping when ready)
; - binary data output
; - disable serial output
; - enable data recorder
CF11101
$W62
; coordinate transformation:
; - radial beam coordinates (2 bits)
; - use pitch/roll (not used for beam coords?)
; - no 3-beam solutions
; - no bin mapping
EX00100
$W62
; Sensor source:
; - manual speed of sound (EC)''
; - manual depth of transducer (ED = 0 [dm])
; - measured heading (EH)
; - measured pitch (EP)
; - measured roll (ER)
; - manual salinity (ES = 35 [psu])
; - measured temperature (ET)
EZ0011101
$W62
; Set transducer depth to zero
ED0000
$W62
; Set salinity to 35ppt
ES35
$W62
;
; - configure no. of bins, length, blank
; number of bins
WN015
$W62
; bin length [cm]
WS0800
$W62
; blank after transmit [cm]
WF0000
$W62
; ambiguity velocity [cm]
WV250
$W62
; amplitude and correlation thresholds for bottom detection
LZ30,220
$W62
; Set ADCP to narrow bandwidth and extend range by 10%
LW1
$W62
; Name data file
RN MLADCP
$W62
;
; Set one ensemble/sec
TE00000100
$W62
; Set one second between pings
TP000100
$W62
; Set LADCP to output Velocity, Correlations, Amplitude, and Percent Good
LD111100000
$W62
; Set one ping per ensemble. Use WP if LADCP option is not enabled.
LP1
$W62
$W62
; keep params as user defaults (across power failures)
CK
$W62
; echo configuration
T?
$W62
W?
$W62
; start Pinging
CS
; Delay 3 seconds
$D3
$p *************************************************************************
$P Please disconnect the ADCP from the computer.
$P *************************************************************************
; Close the log file
$L
```

Recovery:
- Connect cables
- On BBTALK click the terminal window
- press <END> to wake up the ADCP
- In the menu: File --> Recover Recorder --> Select All Files --> OK
- Download selected file following the GUI instructions
- rename the file including the station number, ex: `sta0032_MLADC000.000`
- Put the ADCP to sleep: Type `CZ` and pres <ENTER>

# Post-Processing

## VM-ADCP

Post-processing of this data is necessary to remove biases from the ships movements and misalignment (Angle offset) as well as erroneous backscatter from the seafloor, nets and CTDs, ringing and bubble clouds. We used the CODAS software environment for this purpose.

In the first step, we installed the virtual computer with the CODAS processing software (https://currents.soest.hawaii.edu/docs/adcp_doc/codas_setup/virtual_computer/index.html). Then we copied the folder containing all VM-ADCP files into the shared folder (named after the cruise ID, we used `codas_shared/siarctic2017/enr` here).

Then start the virtual computer, open a terminal inside the folder and enter the following commands. Follow the CODAS GUI instructions when they pop up.

```
mkdir /media/sf_codas_shared/siarctic2017/enr/config
mkdir /media/sf_codas_shared/siarctic2017/enr/uhdas_data

cd /media/sf_codas_shared/siarctic2017/enr/config/
reform_vmdas.py
python vmdas2uhdas.py

proc_starter.py reform_defs.py

cd /media/sf_codas_shared/siarctic2017/enr
adcptree.py os75 --datatype uhdas --cruisename siarctic2017
cd os75
```
Now create a file called `_py.cnt` inside the folder `/media/sf_codas_shared/siarctic2017/enr/os75`
and copy those lines into the file. Change the values where necessary.

```
## python processing

 --yearbase 2017             ## required, for decimal day conversion
                           ##     (year of first data)
 --cruisename siarctic2017     ## *must* match prefix in config dir
 --dbname zzz              ## database name; in adcpdb.  eg. a0918
                           ##
 --datatype uhdas            ## datafile type
 --sonar os75nb              ## specify instrument letters, frequency,
                           ##     (and ping type for ocean surveyors)
 --ens_len  300            ## averages of 300sec duration
                           ##
 --update_gbin             ## required for this kind of processing
 --configtype  python       ## file used in config/ dirctory is python
                           ##

 --max_search_depth 1500   ## try to identify the bottom and eliminate
                           ##    data below the bottom IF topo says
                           ##    the bottom is shallower than 1000m
```

Then run the routine:
```
quick_adcp.py --cntfile q_py.cnt --auto
```
Now check the misalignment amplitude and angle(phase) from the file `cal/watertrk/adcpcal.out` by running:
```
tail -20 cal/watertrk/adcpcal.out
```
To correct for the misalignment, process the data again with the `--rotate_amplitude` and `--rotate_angle   1.80` that you obtained from the previous command.
```
quick_adcp.py --steps2rerun rotate:apply_edit:navsteps:calib --rotate_amplitude 1.0303 --rotate_angle   1.80   --auto
```
Now you are ready to check and eventually edit the ADCP data manually using a great tool called `autoedit`.
```
cd edit
gautoedit.py
```
In this GUI, you can remove bad data that the algorithms could not filter out themselves. Check the entire dataset and press apply for each edit. Once you are done press quit and apply your edits by running this in the terminal:
```
cd ..
quick_adcp.py --steps2rerun apply_edit:navsteps:calib --auto
```
Now you are finally ready to export your data into a netCDF file (found here `contour/os75/`) by running:
```
adcp_nc.py adcpdb  contour/os75  siarctic2017 os75
```

## L-ADCP


# Analysis




After each cruise, the VM-ADCP data was post-processed using the



![](output_n6xlp9.gif)
