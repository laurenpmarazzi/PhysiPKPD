using LightXML
using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics

## first make sure that you have done all the make commands you need to get the right project ready to run

## set path to csv file for saving output
csv_filename = "../Documents/PhysiPKPD/last_dataframe.csv";

## set up which simulations you want to run
nsamps = 5;

dvals = [1,2,5]*exp10.(2:6)'; # drug concentrations
dvals = [0;dvals[:]]; # don't forget to add the control case!

## get all the info from the xml file
xdoc = parse_file("./config/mymodel.xml");
r = root(xdoc);
up = find_element(r, "user_parameters");
d1up = find_element(up, "PKPD_D1_central_increase_on_dose");

overall = find_element(r,"overall")
max_time = find_element(overall,"max_time")
mt = parse(Int64,content(max_time))
csv_datainterval = find_element(up,"csv_data_interval")

## set up output variables
t_vals = 0:parse(Int64,content(csv_datainterval)):mt # time values we'll get cell counts at
DF = DataFrame(time_series=Vector{Int64}[],dose=Float64[]) # the data frame we'll use to plot
export_df = DataFrame(time=[NaN;t_vals]) # the data frame we'll use to export to csv (I know we only want one data frame, but I need to level up my Julia powers first)

mycommand = `./project ./config/mymodel.xml` # the command to run the PhysiPKPD code
cleanup_command = `make data-cleanup` # store the cleanup command here

for dvi = 1:length(dvals)
    set_content(d1up, string(dvals[dvi])) # update the drug concentration
    save_file(xdoc, "./config/mymodel.xml") # save the xml

    for si = 1:nsamps
        run(cleanup_command) # make sure output folder is clean
        run(mycommand) # run code
        df = CSV.File("./output/cell_counts.csv"; header=["time","cell_count"]) |> DataFrame; # store the alive cell counts
        next_col_name = string(dvals[dvi],"-",si) # column names are given by [dose amount]-[sample number]
        export_df[:,next_col_name] = [dvals[dvi];df.cell_count] # append column to the output csv file

        push!(DF,(df.cell_count,dvals[dvi])) # append data to the plotting data frame (including the dose amount)
    end
end

gDF = groupby(DF, :dose) # group data by dose
mean_DF = DataFrame(time_series=Vector{Float64}[],dose=Float64[]) # data frame with mean time series for each dose 
for i = 1:length(gDF)
    push!(mean_DF,(mean(gDF[i].time_series),gDF[i].dose[1]))
end
plot(t_vals,mean_DF.time_series,label = mean_DF.dose') # plot mean time series

CSV.write(csv_filename, export_df;writeheader=false) # write output

## graveyard


# lab_names = string.(mean_DF.dose)
# plot(t_vals,Matrix(mean_DF[:,2:end])',label = string.(mean_DF.dose) )
# plot(t_vals,Matrix(mean_DF[:,2:end])',label = ["100000", "1e6"])
# # CSV.write("../Documents/PhysiPKPD/last_dataframe.csv", DF)

# lab_names2 = ["100000" "1e6"]

# x = 1:10; y = rand(10, 2) # 2 columns means two lines
# plot(t_vals, Matrix(mean_DF[:,2:end])', title = "Two Lines", label = ["100000" "1e6"], lw = 3)
# @df DF scatter(:time,:100.0-1)
# plot(t_vals,Matrix(DF[:,2:end]))