using LightXML
using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics

xdoc = parse_file("./config/mymodel.xml");
r = root(xdoc);
up = find_element(r, "user_parameters");
d1up = find_element(up, "PKPD_D1_central_increase_on_dose");

nsamps = 5;

dvals = [1,2,5]*exp10.(2:6)';
# dvals = [1e5,1e6]

overall = find_element(r,"overall")
max_time = find_element(overall,"max_time")
mt = parse(Int64,content(max_time))
csv_datainterval = find_element(up,"csv_data_interval")
t_vals = 0:parse(Int64,content(csv_datainterval)):mt
DF = DataFrame(time_series=Vector{Int64}[],dose=Float64[])
export_df = DataFrame(time=[NaN;t_vals])

cleanup_command = `make data-cleanup`

for dvi = 1:length(dvals)
    set_content(d1up, string(dvals[dvi]))
    save_file(xdoc, "./config/mymodel.xml")

    for si = 1:nsamps
        run(cleanup_command)
        mycommand = `./project ./config/mymodel.xml`
        run(mycommand)
        df = CSV.File("./output/cell_counts.csv"; header=["time","cell_count"]) |> DataFrame;
        next_col_name = string(dvals[dvi],"-",si)
        export_df[:,next_col_name] = [dvals[dvi];df.cell_count]

        push!(DF,(df.cell_count,dvals[dvi]))
    end
end

gDF = groupby(DF, :dose)
mean_DF = DataFrame(time_series=Vector{Float64}[],dose=Float64[])
for i = 1:length(gDF)
    push!(mean_DF,(mean(gDF[i].time_series),gDF[i].dose[1]))
end
plot(t_vals,mean_DF.time_series,label = mean_DF.dose')

CSV.write("../Documents/PhysiPKPD/last_dataframe.csv", export_df;writeheader=false)
# lab_names = string.(mean_DF.dose)
# plot(t_vals,Matrix(mean_DF[:,2:end])',label = string.(mean_DF.dose) )
# plot(t_vals,Matrix(mean_DF[:,2:end])',label = ["100000", "1e6"])
# # CSV.write("../Documents/PhysiPKPD/last_dataframe.csv", DF)

# lab_names2 = ["100000" "1e6"]

# x = 1:10; y = rand(10, 2) # 2 columns means two lines
# plot(t_vals, Matrix(mean_DF[:,2:end])', title = "Two Lines", label = ["100000" "1e6"], lw = 3)
# @df DF scatter(:time,:100.0-1)
# plot(t_vals,Matrix(DF[:,2:end]))