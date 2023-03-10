image_out_dir="/home/cadel/projects/thesis/UNSWThesis/images"
pres_out_dir="/home/cadel/projects/thesis/UNSWThesis/pres"

sns_style = "darkgrid"
sns_context = "paper"

def save_facet_grid(fg):
    fg.savefig(os.path.join(image_out_dir, "threads_vs_perf.eps"))

assemblies = ["canetoad", "tigersnake", "sandy", "starling"]
excluded = ["canetoad.v3.7b.repeats", "canetoad.v3.7d.junk", "canetoad.v3.7c.quarantine", 
            "canetoad.v3.7a.core", "canetoad.v3.8b.plusmt","sandy.v2.2a.plusmt", "BUSCOMP"]
