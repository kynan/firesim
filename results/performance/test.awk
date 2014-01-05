/Domain/ { sz = $4 }
/MLUP/ { mlup = $3; fmlup = $4 }
END { print sz, mlup, fmlup }
