#!/usr/bin/env python
"""Read an ISC events list in csv format

   This module provides tools to read an 
   list of events from the ISC and turn
   this into a dictionary of event picks
   for a given phase. 

   When run as a script, this ...
"""

import datetime

def _make_datetime(date, time):
    """Date and time are strings"""
    yr, mo, dy = date.split('-', 3)
    hr, mi, se = time.split(':', 3)
    # Sometimes we don't have a decimal
    if (len(se) == 2):
        mse = 0
    else:
        se, mse = se.split('.',2)
    # We assume we have 2 dp of seconds,
    # and convert the value (0 - 99) to
    # microseconds
    dati = datetime.datetime(int(yr), int(mo), 
                             int(dy), int(hr), 
                             int(mi), int(se),
                             int(mse)*10000 )
    return dati

def read_picks(filename, phaselist):

    fh = open(filename, 'r')
    # lists, for the types of pick
    all_picks = {}
    for phase in phaselist:
        all_picks[phase] = {}

    for line in fh:
        # New dictionary for this pick
        thispick = {}
        words = line.split(',')
        thispick['eventid'] = words[0].strip()
        thispick['reporter'] = words[1].strip()
        thispick['station'] = words[2].strip()
        thispick['station_lat'] = float(words[3].strip())
        thispick['station_lon'] = float(words[4].strip())
        thispick['station_elev'] = float(words[5].strip())
        thispick['epicentral_distance'] = float(words[7].strip())
        thispick['backazimuth'] = float(words[8].strip())
        thispick['phase'] = words[9].strip()
        thispick['pick_datetime'] = _make_datetime(words[11].strip(), words[12].strip())
        thispick['event_datetime'] = _make_datetime(words[18].strip(), words[19].strip())
        thispick['event_lat'] = float(words[20].strip())
        thispick['event_lon'] = float(words[21].strip())
        thispick['event_depth'] = float(words[22].strip())
        # Sometimes the same reporter reports the same phase at the same station
        # more than once. We just use the most recent. 
        # We return the results in a dict (of dicts, of dicts)
        pick_key = thispick['eventid']+thispick['station']+thispick['reporter']
        for phase in phaselist:
            if thispick['phase'] == phase:
                all_picks[phase][pick_key] = thispick             
    fh.close()
    
    return all_picks

def pair_picks(all_picks, phase1, phase2):

    pick_pairs = {}
    count = 0
    # Find the P's (phase1)  given the PcP's (phase 2)(if we have them)
    for event_station in all_picks[phase1]:
        if event_station in all_picks[phase2]:
            pick1 = all_picks[phase1][event_station]
            pick2 = all_picks[phase2][event_station]
            thispick = {}
            thispick['event_lat'] = pick1['event_lat']
            thispick['event_lon'] = pick1['event_lon']
            thispick['event_depth'] = pick1['event_depth']
            thispick['event_datetime'] = pick1['event_datetime']
            thispick['station'] = pick1['station']
            thispick['reporter'] = pick1['reporter']
            thispick['eventid'] = pick1['eventid']
            thispick['station_lat'] = pick1['station_lat']
            thispick['station_lon'] = pick1['station_lon']
            thispick['station_elev'] = pick1['station_elev']
            thispick['epicentral_distance'] = pick1['epicentral_distance']
            thispick['backazimuth'] = pick1['backazimuth']

            thispick[phase1+'_datetime'] = pick1['pick_datetime']
            thispick[phase2+'_datetime'] = pick2['pick_datetime']
            pick_pairs[event_station] = thispick

    return pick_pairs
