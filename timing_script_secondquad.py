"""
To be called like this:

python -m cProfile timing_script_secondquad.py

"""

print "Good morning."

from astrodendro_analysis.second_quadrant_cloud_extraction import second_quad_dendrogram

print "I have imported the second_quad_dendrogram function."
print "Now I am going to call it:"

d2, catalog2, header2, metadata2 = second_quad_dendrogram()

print "I have finished."