#  =======================================================================
#           R-ScaLAPACK version 0.3.x:  R interface to ScaLAPACK
#              Oak Ridge National Laboratory, Oak Ridge TN.
#      Authors: David Bauer, Nagiza. F. Samatova, Srikanth Yoginath
#     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
#                 Computer Science and Mathematics Division
#             Oak Ridge National Laboratory, Oak Ridge TN 37831 
#                   (C) 2004 All Rights Reserved
#
#                              NOTICE
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby granted
# provided that the above copyright notice appear in all copies and
# that both the copyright notice and this permission notice appear in
# supporting documentation.
#
# Neither the Oak Ridge National Laboratory nor the Authors make any
# representations about the suitability of this software for any
# purpose.  This software is provided ``as is'' without express or
# implied warranty.
#
# R-ScaLAPACK (http://www.aspect-sdm.org/R-ScaLAPACK) was funded
# as part of the Scientific Data Management Center
# (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
# Scientific Discovery through Advanced Computing (DOE SciDAC) program
# (http://www.scidac.org ). 
# ========================================================================

PA.exit <- function (){
	.Call("PA_Exit", PACKAGE="RScaLAPACK")
}

PA.exec <- function (tempPackage, tempScript, inputVector){

	# The task of this function is to execute the specified 
	# R script in parallel and wait for the result.
	# On obtaining the result back from the spawned processes,
	# it should return it to the called function.


	scriptLocn <- system.file("exec",tempScript,package=tempPackage)

	x<- .Call("PA_Exec", as.character(scriptLocn), inputVector,
		PACKAGE="RScaLAPACK")


	return(x);

}
