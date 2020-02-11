%ROWINDEX = ROWFIND(VALUES, LOOKUPVALUES)
%
%Finds rows of LOOKUPVALUES that are equal to rows in VALUES and return the index
%of the rows.  If more than one row from LOOKUPVALUES matches, it returns
%the index of the first matching row only. If no row matches, it returns a
%zero for that row.  The number of columns in VALUES and LOOKUPVALUES must
%be equal.
%
%