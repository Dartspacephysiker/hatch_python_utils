# 2019/07/03
import numpy as np
from datetime import datetime

# Fra Kjellmar-epost 20190702
# The header lines contain the following information:
#   Date
#   Start_time (Universal Time)
#   End_time (Universal Time)
#   Antenna_azimuth (degrees east of geographic north)
#   Antenna_elevation (degrees above horizon)
#   Geographic latitude of radar site (degrees)
#   Geographic longitude of radar site
#   Altitude of radar site (km above sea level)
#   Number of data lines to follow (in the case shown above: 26 data lines. This number will vary with time, and should be helpful when you write software to load the data).
#
# Each of the data lines contains the following information:
#   Height of observations (km)
#   Geographic Latitude (deg)
#   Geographic Longitude (deg)
#   Log10 of the electron density (m^-3)
#   Electron temperature (Kelvins)
#   Ion temperature (Kelvins)
#   Ion velocity (m/s, positive values are upflow, i.e. away from the radar)
#   Error of electron density (measured in percent)
#   Error of electron temperature (measured in percent)
#   Error of ion temperature (measured in percent)
#   Error of ion velocity (measured in percent)
#


def Txt_Reader(dataDir, dataFile,
               verbose=False):

    recDType = np.dtype([('recno', int),
                         ('startT', datetime),
                         ('stopT', datetime),
                         ('az', float),
                         ('el', float),
                         ('radarLat', float),
                         ('radarLon', float),
                         ('radarAlt', float),
                         ('obsAlt', float),
                         ('obsLat', float),
                         ('obsLon', float),
                         ('Nel_log10', float),
                         ('Tel', float),
                         ('Tion', float),
                         ('Vion', float),
                         ('Nel_err_pct', float),
                         ('Tel_err_pct', float),
                         ('Tion_err_pct', float),
                         ('Vion_err_pct', float)])

    with open(dataDir+dataFile) as f:
        lines = f.readlines()

        nLines = len(lines)

        records = np.recarray((nLines,), dtype=recDType)

        iLine = 0
        nRec = 0
        totLines = 0
        while iLine < nLines:

            # Read line
            while lines[iLine].isspace():

                iLine += 1

                if iLine == nLines:
                    break

            if iLine == nLines:
                break

            lineL = [thing for thing in lines[iLine].split(
                ' ') if ((thing != '') and (thing != '\n'))]

            # Start_time and End_time (Universal Time)
            startT, stopT = lineL[1].split('-')
            startT = datetime.strptime(
                lineL[0]+' '+startT, "%d-%b-%Y %H:%M:%S")
            stopT = datetime.strptime(lineL[0]+' '+stopT, "%d-%b-%Y %H:%M:%S")

            # Antenna_azimuth (degrees east of geographic north)
            az = float(lineL[2])
            el = float(lineL[3])  # Antenna_elevation (degrees above horizon)
            # Geographic latitude of radar site (degrees)
            radarLat = float(lineL[4])
            radarLon = float(lineL[5])  # Geographic longitude of radar site
            # Altitude of radar site (km above sea level)
            radarAlt = float(lineL[6])
            nDataLines = int(lineL[7])  # Number of data lines to follow

            iLine += 1

            if verbose:
                print("Reading {:d} lines for record #{:d} ...".format(
                    nDataLines, nRec))
            iDataLine = 0
            while iDataLine < nDataLines:
                dataLine = [thing for thing in lines[iLine +
                                                     iDataLine].split(' ') if ((thing != '') and (thing != '\n'))]

                obsAlt = float(dataLine[0])  # Height of observations (km)
                obsLat = float(dataLine[1])  # Geographic Latitude (deg)
                obsLon = float(dataLine[2])  # Geographic Longitude (deg)
                # Log10 of the electron density (m^-3)
                Nel_log10 = float(dataLine[3])
                Tel = float(dataLine[4])  # Electron temperature (Kelvins)
                Tion = float(dataLine[5])  # Ion temperature (Kelvins)
                # Ion velocity (m/s, positive values are upflow, i.e. away from the radar)
                Vion = float(dataLine[6])
                # Error of electron density (measured in percent)
                Nel_err_pct = float(dataLine[7])
                # Error of electron temperature (measured in percent)
                Tel_err_pct = float(dataLine[8])
                # Error of ion temperature (measured in percent)
                Tion_err_pct = float(dataLine[9])
                # Error of ion velocity (measured in percent)
                Vion_err_pct = float(dataLine[10])

                records[totLines] = (nRec, startT, stopT, az, el, radarLat, radarLon, radarAlt,
                                     obsAlt, obsLat, obsLon,
                                     Nel_log10, Tel, Tion, Vion,
                                     Nel_err_pct, Tel_err_pct, Tion_err_pct, Vion_err_pct)

                totLines += 1

                iDataLine += 1

                # break

            nRec += 1
            iLine += nDataLines

    return records[0:totLines]
