#! /bin/sh

BASE_PATH=`echo $0 | sed "s/\(.*\/\)[^/]*$/\1/"`

if test -e "$BASE_PATH/../libs/RScaLAPACK.so"; then
	export LD_LIBRARY_PATH="$BASE_PATH/../libs:$LD_LIBRARY_PATH"
else
	LOCATE_LIB_LONG=`locate /RScaLAPACK.so | sed -n "s/\(.*\/\)[^/]*$/\1/p"`
	LOCATE_LIB_LONG=`echo $LOCATE_LIB_LONG | sed "s/[ ]/:/g"`
	export LD_LIBRARY_PATH=$LOCATE_LIB_LONG:$LD_LIBRARY_PATH
fi

if test -e "$BASE_PATH/CRDriver"; then
	$BASE_PATH/CRDriver
else
	EXEC_LONG=`which CRDriver`
	EXEC_SHORT=`echo $EXEC_LONG | sed -n "s/.*\/\([^/]*\)$/\1/p"`

	if test -n "$EXEC_SHORT" && test "CRDriver" == $EXEC_SHORT; then
		CRDriver
	else
		LOCATE_DRIVERS=`locate /exec/CRDriver`
		NUM_DRIVERS=`echo $LOCATE_DRIVERS | wc -l`
		if test $NUM_DRIVERS == "1"; then
			$LOCATE_DRIVERS
		fi
	fi

# TODO More tests
fi	

