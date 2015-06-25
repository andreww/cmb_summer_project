#!/usr/bin/env python

"""Download .csv file from the ISC.

   This module provides tools to Download
   .csv formatted files of phase arrival
   times from the ISC catalogue. It does 
   this by formatting a URL and sending a 
   download request for the specified data
   to the ISC
"""
import urllib

def get_csv(url):
    urllib.urlretrieve(url,'isc_data.csv')
    return

def url_build(url_dict):
    """Builds a url from a dictionary containing the search parameters"""
    isc_url = str(url_dict['address']+url_dict['out_format']+'&'+
                  url_dict['request_type']+'&'+url_dict['station_region']+'&'
                  +url_dict['arrivals_limits']+'&'+url_dict['phaselist']+'&'+
                  url_dict['event_region']+'&'+url_dict['start_time']+'&'+
                  url_dict['end_time'])
    return isc_url
# Create a quick test dictionary for building url
url_params = {}
url_params['address'] = 'http://isc-mirror.iris.washington.edu/cgi-bin/web-db-v4?'
url_params['out_format'] = 'out_format=CSV'
url_params['request_type'] = 'request=STNARRIVALS'
url_params['arrivals_limits'] = 'ttime=on&iscreview=on'
url_params['station_region'] = 'sta_list&stnsearch=GLOBAL'
url_params['event_region'] = 'searchshape=GLOBAL'
url_params['phaselist'] = 'phaselist=P,PcP'
url_params['start_time'] = 'start_year=2011&start_month=8&start_day=01&start_time=00:00:00'
url_params['end_time'] = 'end_year=2012&end_month=1&end_day=01&end_time=00:00:00'

url = url_build(url_params)
print url
get_csv(url)