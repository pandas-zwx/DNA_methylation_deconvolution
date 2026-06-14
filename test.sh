mkdir -p /tmp/bedtools_map_test
cd /tmp/bedtools_map_test

cat > CpG_clean.bed <<'EOF'
chr1	100	102	CpG_1
chr1	200	202	CpG_2
chr1	300	302	CpG_3
EOF

cat > bed_file.bed <<'EOF'
chr1	100	102	0.5
chr1	300	302	0.8
EOF

echo '=== CpG_clean.bed ==='
cat CpG_clean.bed
echo

echo '=== bed_file.bed ==='
cat bed_file.bed
echo

echo '=== test: without -null ==='
bedtools map \
  -a CpG_clean.bed \
  -b bed_file.bed \
  -c 4 \
  -o distinct \
  -f 1.0 \
  -F 1.0
echo

echo '=== test: with -null -1 ==='
bedtools map \
  -a CpG_clean.bed \
  -b bed_file.bed \
  -c 4 \
  -o distinct \
  -null -1 \
  -f 1.0 \
  -F 1.0
echo

echo '=== final pipeline ==='
bedtools map \
  -a CpG_clean.bed \
  -b bed_file.bed \
  -c 4 \
  -o distinct \
  -null -1 \
  -f 1.0 \
  -F 1.0 \
| awk '$5 >= 0 {print $1, $2, $3, $5}' OFS='\t'

