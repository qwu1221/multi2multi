# A multivariate to multivariate approach for voxel-wise genome-wide association analysis

"bipar_greedy.m" is from Algorithm 1 to maximize the objective function (4)

"greedy_lik.m" is to select the tunning parameter via likelihood function using Algorithm A1

"Refine_spatial.m" corresponds to Algorithm A2 that extract IGDB with spatial constraint. Refine_spatial.m needs Dijkstra Algorithms by "Dimas Aryo (2021). Dijkstra Algorithm (https://www.mathworks.com/matlabcentral/fileexchange/36140-dijkstra-algorithm), MATLAB Central File Exchange. Retrieved November 2, 2021."

"demo_multiple_IGDB.m" is a simulation example when multiple IGDBs appear. Its correponding results are included as html file. 

"demp_spatial.m" illustrats the implementation of "Refine_spatial.m" using simulated data. The matlab published results are included in html.
