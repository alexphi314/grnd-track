import requests
import os

QUERY_URL = 'https://www.space-track.org/basicspacedata/query/class/tle_latest/NORAD_CAT_ID/{}/orderby/ORDINAL asc/metadata/false'
LOGIN_URL = 'https://www.space-track.org/ajaxauth/login'

def get_tle(cat_id):
    """
    Given cat id, return the most recent TLE for that object

    :param cat_id: String or double
    :return: TLE string
    """

    username = os.getenv('SATCAT_USER')
    password = os.getenv('SATCAT_PASSWORD')
    if username is None or password is None:
        raise Exception('SATCAT_USER and SATCAT_PASSWORD must be defined credentials in environment!')

    payload = {'identity': username, 'password': password}

    r = requests.post(LOGIN_URL, data=payload)
    if r.status_code != 200:
        raise Exception('Issue connecting to space-track login. Status code: {}'.format(r.status_code))

    del password
    del payload

    sat_cookies = r.cookies
    r = requests.get(QUERY_URL.format(cat_id), cookies=sat_cookies).json()

    line1 = r[0]['TLE_LINE1']
    line2 = r[0]['TLE_LINE2']

    return line1, line2