#!/usr/bin/env python3
from __future__ import print_function
import requests
from bs4 import BeautifulSoup
import re
import os
import urllib3
urllib3.disable_warnings()


def getData():
    url = 'https://strike.scec.org/cvws/cgi-bin/seas.cgi'
    stations = ['fltst_dp000', 'fltst_dp025', 'fltst_dp050', 'fltst_dp075', 
                'fltst_dp100', 'fltst_dp125', 'fltst_dp150', 'fltst_dp175', 
                'fltst_dp200', 'fltst_dp250', 'fltst_dp300', 'fltst_dp350']
    with requests.Session() as session:
        res = session.get(url, verify=False)
        soup = BeautifulSoup(res.content, "html.parser")
        res = session.post(url, data={'o': '1005', 'G0012': 'Go -->'})
        res = session.post(url, data={'o': '1005', 'G1045bp1-qd': ' Select '})
        # print(res.content)
        soup = BeautifulSoup(res.content, "html.parser")
        groups = soup.find_all(
            'input', {'type': 'submit', 'value': ' Select '})
        for group in groups:
            author = group.get('name')[5::]
            print('* ' + author)
            res = session.post(
                url, data={'m': 'bp1-qd', 'o': '1005', group.get('name'): ' Select '})
            soup = BeautifulSoup(res.content, "html.parser")
            files = soup.find_all(
                'input', {'type': 'submit', 'value': 'Raw Data'})
            mus = soup.find(
                'input', {'type': 'hidden', 'name': 'mus'}).get('value')
            for file in files:
                file_station = file.get('name')[5::]
                if file_station in stations:
                    print('\t' + file_station)
                    res = session.post(
                        url, data={'m': 'bp1-qd', 'o': '1005', 'mus': mus, file.get('name'): 'Raw Data'})
                    soup = BeautifulSoup(res.content, "html.parser")
                    data = soup.find('pre')
                    if not os.path.exists(file_station):
                        os.makedirs(file_station)
                    f = open('%s/%s.txt' % (file_station, author), 'w')
                    f.writelines(data.contents)
                    f.close()


if __name__ == '__main__':
    getData()
