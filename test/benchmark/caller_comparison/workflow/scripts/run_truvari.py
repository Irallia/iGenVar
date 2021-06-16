import subprocess

subprocess.run(["rm", "-rf", snakemake.params["out_dir"]])
if snakemake.params["vcf"].endswith(".gt.seq"):
	process = subprocess.Popen(["truvari", "-f", snakemake.input["genome"], "-b", snakemake.input["truth_vcf"], 
								"-c", snakemake.input["calls"], "-o", snakemake.params["out_dir"], "--passonly",
								"--includebed", snakemake.input["truth_bed"], "--sizemax", "1000000", "-r", "1000", 
								"--pctsim", "0.7", "--pctsize", "0.7", "--pctovl", "0", "--gtcomp"], stderr=subprocess.PIPE)
elif snakemake.params["vcf"].endswith(".gt"):
	process = subprocess.Popen(["truvari", "-f", snakemake.input["genome"], "-b", snakemake.input["truth_vcf"], 
								"-c", snakemake.input["calls"], "-o", snakemake.params["out_dir"], "--passonly",
								"--includebed", snakemake.input["truth_bed"], "--sizemax", "1000000", "-r", "1000", 
								"--pctsim", "0", "--pctsize", "0.7", "--pctovl", "0", "--gtcomp"], stderr=subprocess.PIPE)
elif snakemake.params["vcf"].endswith(".seq"):
	process = subprocess.Popen(["truvari", "-f", snakemake.input["genome"], "-b", snakemake.input["truth_vcf"], 
								"-c", snakemake.input["calls"], "-o", snakemake.params["out_dir"], "--passonly",
								"--includebed", snakemake.input["truth_bed"], "--sizemax", "1000000", "-r", "1000", 
								"--pctsim", "0.7", "--pctsize", "0.7", "--pctovl", "0"], stderr=subprocess.PIPE)
else:
	process = subprocess.Popen(["truvari", "-f", snakemake.input["genome"], "-b", snakemake.input["truth_vcf"], 
								"-c", snakemake.input["calls"], "-o", snakemake.params["out_dir"], "--passonly",
								"--includebed", snakemake.input["truth_bed"], "--sizemax", "1000000", "-r", "1000", 
								"--pctsim", "0", "--pctsize", "0.7", "--pctovl", "0"], stderr=subprocess.PIPE)

# wait for the process to terminate
out, err = process.communicate()

with open(snakemake.log["log"], 'w') as log_file:
	print(err, file=log_file)
if "[INFO] Finished" in str(err):
	print("Truvari finished successfully")
elif "Specific error: \"fetch requires an index\"" in str(err):
	print("Truvari failed on empty VCF file")
	with open(snakemake.output["summary"], 'w') as summary:
		print("{", file=summary)
		print("\t\"precision\": 0,", file=summary)
		print("\t\"recall\": 0", file=summary)
		print("}", file=summary)
else:
	raise Exception(err)