import pandas as pd
from geojson import MultiLineString
import scipy.io

datafiles = ['Zone1_Soln_Route.mat', 'Zone4_Soln_Route.mat', 'Zone5_Soln_Route.mat']

points_arr = []
for dfile in datafiles:
    df = scipy.io.loadmat(dfile)['route_long_lat']
    df_clip = pd.DataFrame([df[i] for i in range(len(df)) if df[i][0] != 0])
    points = [(df_clip.iloc[i,0], df_clip.iloc[i,1]) for i in range(df_clip.shape[0])]
    points_arr.append(points)

route_str = MultiLineString(points_arr)

f = open("routes.json", "w")
f.write(str(route_str))
f.close()
