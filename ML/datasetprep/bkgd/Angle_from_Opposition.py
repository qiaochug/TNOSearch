
import argparse
import ephem
import numpy as np
import pandas as pd

out = pd.DataFrame()

def add_angle(inputcsv):
    global out
    
    df = pd.read_csv(inputcsv)
    df.insert(4,'ang_from_opp', pd.Series(np.zeros(len(df['SNOBJID'])), df.index))
    for index, row in df.iterrows():
        lon,lat,angle = angle_from_opp(row['RA'], row['DEC'], row['MJD'])
        df.loc[index, 'ang_from_opp']= angle  
    out = out.append(df, ignore_index = True)

# This function needs the ephem python package
def angle_from_opp(ra,dec,date):
    '''
    ra, dec are ephem.Angle() objects
    date is an ephem.date() object
    output value from -pi to pi
    '''
    ra = ephem.hours(ra)
    dec = ephem.degrees(dec)
    date = ephem.Date(date - 15020)
    
    lon, lat = ephem.Ecliptic(ephem.Equatorial(ra,dec)).get()
    sun = ephem.Sun()
    sun.compute(date)
    lon_s, lat_s = ephem.Ecliptic(ephem.Equatorial(sun.ra, sun.dec)).get()
    delta_lon = ephem.degrees(lon_s - lon).znorm
    
    ang_from_opp = np.pi - delta_lon
    if ang_from_opp < - np.pi:
        ang_from_opp = ang_from_opp + 2*np.pi
    if ang_from_opp > np.pi:
        ang_from_opp = ang_from_opp - 2*np.pi
    
    if ang_from_opp > np.pi or ang_from_opp < -np.pi:
        print("Abnormal value!!")
    
    return lon, lat, ang_from_opp

def main():
    """
    takes in an output path csv
    takes in the format of the input csv
    transform each input csv by calculating opp for each entry and adding them in the dataframe,
    then append the result csv to the end of the outputcsv
    """
    global out
    parser = argparse.ArgumentParser()
    parser.add_argument('output', nargs = 1, default = None,
                         help = 'output csv file with unique ids')
    parser.add_argument('csvfileinput', nargs = 1, default = None,
                         help = 'input csvfils')
    args = parser.parse_args()
    for num in np.arange(1,41):
        add_angle(args.csvfileinput[0]+ str(num) + '.csv')
    out.to_csv(args.output[0])
    
if __name__ == '__main__':
    main()