docker run -it \
  -v ~/Desktop/test_data:/in_data \
  -v ~/Desktop/test_data:/out_data
  sushi/rm_data_module:latest \
  -d /in_data/normalized_counts.csv \
  -o /out_data/outfile.csv \
  -c /in_data/remove_criteria.csv
