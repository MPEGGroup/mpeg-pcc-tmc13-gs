#!/bin/bash
#
# Generate a configuration tree in $PWD from YAML files in the same
# directory.

set -e
shopt -s nullglob

script_dir="$(cd "$(dirname "$0")" && pwd)"
geom="octree"
attr="predlift"
src_cfg_dir=""
mode="pointcloud"
gs_set=""
frame_id=""

while (( $# )); do
	case $1 in
	--octree) geom="octree" ;;
	--predgeom) geom="predgeom" ;;
	--trisoup) geom="trisoup" ;;
	--raht) attr="raht" ;;
	--pred-lift) attr="predlift" ;;
	--hybrid) attr="hybrid" ;;
	--all) all=1 ;;
	--pointcloud) mode="pointcloud" ;;
	--gaussiansplat) mode="gaussiansplat" ;;
	# Carry a per-frame id as a uint16 reflectance attribute (ConcatGPCC).
	--frame-id) frame_id=1 ;;
	--1F)
		[[ -z "$gs_set" || "$gs_set" == "1F" ]] || { echo "error: cannot combine --1F and --NF" >&2; exit 1; }
		gs_set="1F"
		;;
	--NF)
		[[ -z "$gs_set" || "$gs_set" == "NF" ]] || { echo "error: cannot combine --1F and --NF" >&2; exit 1; }
		gs_set="NF"
		;;
	--cfgdir=*) src_cfg_dir="${1#--cfgdir=}/" ;;
	--) shift; break ;;
	--help|*)
		echo -e "usage:\n $0\n" \
			"    [[--octree|--predgeom|--trisoup] [--raht|--pred-lift|--hybrid] | --all]\n" \
			"    [--pointcloud|--gaussiansplat [--1F|--NF] [--frame-id]]\n" \
			"    [--cfgdir=<dir>]"
		exit 1
	esac
	shift;
done

extra_args=("$@")

##
# NB: it is important that the configs in each config set are
# capable of being merged together by gen-cfg.pl.  Ie, no two
# configs may have different definitions of one category.
cfg_octree_predlift=(
	octree-liftt-ctc-lossless-geom-lossy-attrs.yaml
	octree-liftt-ctc-lossy-geom-lossy-attrs.yaml
	octree-predt-ctc-lossless-geom-lossless-attrs.yaml
	octree-predt-ctc-lossless-geom-nearlossless-attrs.yaml
)

cfg_octree_raht=(
	octree-raht-ctc-lossless-geom-lossy-attrs.yaml
	octree-raht-ctc-lossy-geom-lossy-attrs.yaml
)

cfg_octree_hybrid=(
	octree-hybrid-ctc-lossless-geom-lossy-attrs.yaml
)

cfg_predgeom_predlift=(
	"${cfg_octree_predlift[@]}"
	cfg-predgeom.yaml
)

cfg_predgeom_raht=(
	"${cfg_octree_raht[@]}"
	cfg-predgeom.yaml
)

cfg_predgeom_hybrid=(
	"${cfg_octree_hybrid[@]}"
	cfg-predgeom.yaml
)

cfg_trisoup_predlift=(
	trisoup-liftt-ctc-lossy-geom-lossy-attrs.yaml
)

cfg_trisoup_raht=(
	trisoup-raht-ctc-lossy-geom-lossy-attrs.yaml
)

# ConcatGPCC: with --frame-id, use the config variants that add a uint16
# reflectance attribute carrying the per-frame id.  Only the two CTC GS
# configs (octree-hybrid C1, octree-predlift CW) are needed for that pipeline.
if [[ -n "$frame_id" ]]; then
	cfg_octree_hybrid=(
		octree-hybrid-ctc-lossless-geom-lossy-attrs-frameid.yaml
	)
	cfg_octree_predlift=(
		octree-predt-ctc-lossless-geom-lossless-attrs-frameid.yaml
	)
fi

do_one_cfgset() {
	local geom=$1
	local attr=$2

	outdir="$geom-$attr/"
	mkdir -p "$outdir"

	local cfgset_name
	local cfgset_len
	cfgset_name="cfg_${geom}_${attr}"
	cfgset="${cfgset_name}[@]"

	eval "cfgset_len=\${#${cfgset_name}[@]}"
	if (( cfgset_len == 0 )); then
		echo "No config set for mode=$mode geom=$geom attr=$attr" >&2
		return 1
	fi

	for f in ${!cfgset}
	do
		echo "${src_cfg_dir}$f -> $outdir" ...
	done

	sequence_yaml_args=()
	mode_extra_args=()
	if [[ "$mode" == "pointcloud" ]]; then
		sequence_yaml_args=(
			"${src_cfg_dir}sequences-cat1.yaml"
			"${src_cfg_dir}sequences-cat3.yaml"
		)

		# Ensure catGS sequences does not get generated in pointcloud mode.
		if [[ -f "${src_cfg_dir}sequences-catGS.yaml" ]]; then
			gs_sequences=($(grep -E '^[[:space:]]+[a-zA-Z0-9_-]+:' "${src_cfg_dir}sequences-catGS.yaml" \
				| sed -E 's/^[[:space:]]*([a-zA-Z0-9_-]+):.*/\1/'))
			if (( ${#gs_sequences[@]} )); then
				# gen-cfg.pl uses --exlude-seqs (sic)
				mode_extra_args+=(--exlude-seqs="$(IFS=:; echo "${gs_sequences[*]}")")
			fi
		fi
	else
		 sequence_yaml_args=("${src_cfg_dir}sequences-catGS.yaml")

		# Extract sequence names dynamically from YAML
		gaussian_sequences=($(grep -E '^[[:space:]]+[a-zA-Z0-9_-]+:' "${src_cfg_dir}sequences-catGS.yaml" \
			| sed -E 's/^[[:space:]]*([a-zA-Z0-9_-]+):.*/\1/'))

		# Default GS set is 1F (all sequences in sequences-catGS.yaml)
		: "${gs_set:=1F}"

		if [[ "$gs_set" == "NF" ]]; then
			filtered_sequences=()
			for s in "${gaussian_sequences[@]}"; do
				case "$s" in
					bartender_*|cinema_*|breakfast_*|manwithfruit_*|breakdance_*|trio*|makeup*|musical*)
						filtered_sequences+=("$s")
						;;
				esac
			done
			gaussian_sequences=("${filtered_sequences[@]}")
		fi

		if (( ${#gaussian_sequences[@]} == 0 )); then
			echo "error: no gaussiansplat sequences selected (gs_set=$gs_set)" >&2
			exit 1
		fi

		# Pass them to gen-cfg.pl
		mode_extra_args=(
			--only-seqs="$(IFS=:; echo "${gaussian_sequences[*]}")"
		)
	fi

	# NB: specifying extra_args at the end does not affect option
	# processing since gen-cfg.pl is flexible in argument positions
	$script_dir/gen-cfg.pl \
		--prefix="$outdir" --no-skip-sequences-without-src \
		"${!cfgset/#/${src_cfg_dir}}" \
		"${sequence_yaml_args[@]}" \
		"${mode_extra_args[@]}" \
		"${extra_args[@]}"

	# In gaussiansplat mode, prune any legacy outputs that may already exist.
	# (gen-cfg.pl won't touch mtimes if content is unchanged, but it also won't
	# delete obsolete variants / files.)
	if [[ "$mode" == "gaussiansplat" ]]; then
		find "$outdir" -type d -name "r05" -prune -exec rm -rf {} + 2>/dev/null || true
		find "$outdir" -type d -name "r06" -prune -exec rm -rf {} + 2>/dev/null || true
		find "$outdir" -type f -name "pcerror.cfg" -exec rm -f {} + 2>/dev/null || true
	fi

	rm -f "$outdir/config-merged.yaml"
}

if [[ "$mode" == "pointcloud" ]]; then
	if [[ "$all" != "1" ]]
	then
		do_one_cfgset "$geom" "$attr"
	else
		do_one_cfgset "octree" "predlift"
		do_one_cfgset "octree" "raht"
		do_one_cfgset "predgeom" "predlift"
		do_one_cfgset "predgeom" "raht"
		do_one_cfgset "trisoup" "predlift"
		do_one_cfgset "trisoup" "raht"
	fi
else
	if [[ "$all" != "1" ]]
	then
		do_one_cfgset "$geom" "$attr"
	else
		do_one_cfgset "octree" "predlift"
		do_one_cfgset "octree" "raht"
		do_one_cfgset "octree" "hybrid"
		do_one_cfgset "predgeom" "predlift"
		do_one_cfgset "predgeom" "raht"
		do_one_cfgset "predgeom" "hybrid"
	fi
fi
