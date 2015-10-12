#!/bin/bash
# Make sure that the csh version and Fortran 90 module of tomo_predict
# give the same answer for a given 3D Earth model.

# We need some programs to carry on
for p in gfortran taup_path; do
	command -v "$p" >/dev/null 2>&1 || { echo "Need program $p to run test" >&2; exit 1; }
done

# source-receiver geometry and phase name
evlo=0
evla=0
stlo=90
stla=0
phase=PcP

# Input to the csh version: stla stlo evla evlo evdp phase
input="$stla $stlo $evla $evlo $phase"

# 1D and 3D model files
model1d=ak135.1D_vp
model3d=vdh3D_1999

# taup_path output
path=taup_path.gmt

# Run csh script, which creates taup_path.gmt
./c.tomo_predict_vdh "$phase" >/dev/null &&
[ -f "$path" ] || { echo "Some problem running script or creating \"$path\"" >&2; exit 1; }
# Get dt from output
dt_ejg=$(awk 'NR==2{print $1}' vdh1999.$phase)

# Run precompiled version; capture output, which is dt
dt_ajn=$(
	cat <<-END | ./tomo_predict_vdh | awk '{print $1}'
	$model1d
	$model3d
	$path
	0 0
	END
)

# Compare values
echo  "Original version    Rewritten version"
printf "%16.3f    %17.3f\n" "$dt_ejg" "$dt_ajn"
