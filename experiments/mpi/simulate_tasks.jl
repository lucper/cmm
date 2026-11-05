#!/usr/bin/env julia
#
# Simulates the search on several nodes from the times of its tasks (wildcard configurations)
# measured on 1 node, to see how much of the slowdown on more nodes the uneven tasks explain. As
# main_algo does, the tasks are cut into one block of consecutive tasks per process, of sizes
# differing by at most 1, and the threads of each process take the next task of its block as soon
# as they are free (schedule(dynamic)): each thread takes one at the start, then the next free
# thread the next one. A run takes as long as its last thread, plus a fixed overhead. Each trial
# shuffles the tasks at random; if the input gives the place of each task in the order of the code
# (confN), that order is simulated too.
#
# A run spends its time in four parts, which the output separates (and compares with tasks that
# all take their mean, which shows what the number of tasks per thread alone costs):
# - the tasks, as their sizes: the time of each task divided by the speed of the core that ran it,
#   and the first task of each thread at most the median of the others;
# - the speed of each core: a few cores of a run are 10-40% slower than the others through the
#   whole run (on different cores in each run, so the slowness belongs to the core, not to the
#   task). The speed of a thread is the median of its tasks after the first against the median of
#   all of those (1 with fewer than 3 of them), and a task on it takes its size times that speed.
#   Dynamic scheduling gives a slow core fewer tasks, which tasks that were slow in themselves would
#   not get (see --cores);
# - the warm-up of each thread: the time its first task took beyond that median, which every thread
#   pays once, on any number of nodes;
# - the overhead of each process, outside its tasks, such as building its index and searchers,
#   which no number of nodes shortens (see --overhead).
# Each trial gives the measured speeds and warm-ups, of the same threads, to the threads of each
# node in a random order; the order of the code gives each thread of a process those of the thread
# that ran the same task on 1 node.
#
# From experiments/mpi, with the times of each layout on 1 node and the table of the layouts:
#   julia -t auto simulate_tasks.jl 010-measure/human_ii_timers/*.n1.*.trial1.tsv \
#       --results 020-concat/human_ii.layouts.l8d5.x2.k10000.tsv --plot 030-simulate/human_ii.l8d5
# Several files of times are simulated in turn, each with its own figures (prefix-<layout>).
#
# The input has one task per line, with its time in seconds in the last column, or in the column
# given with --column (from 1). As main_algo writes them with the timers, the lines are
# p<process>t<thread>, conf<place in the order of the code>, then optionally the node, the core
# (cpu17, or cpu17to18 if the thread moved) and the start (start<seconds since the epoch)
# of the task, and its time in seconds; the threads and places are optional, and without the
# threads there is no warm-up. Other fields are ignored. Empty lines and lines starting with # are skipped,
# and so is a first line whose time is not a number (a header).
#
# The layout (processes x threads per node) is --layout, or else that of the processes and threads
# of the input, or else that of the name of the file (.n1.8x16.).
#
# --results reads the times measured on the cluster from a table of 020-concat: the lines of the
# layout, with the slowest search of a process (column 21, in ms), averaged over the trials.
# --measured gives them instead as nodes:seconds,... With one on the first number of nodes, the
# measured speedup is compared with the hypothetical one, and the output says how much of the
# lost speedup (perfect minus measured) the hypothetical one explains. The numbers of nodes are
# --nodes, or else those measured, or else 1,2,4,8.
#
# --overhead is the overhead in seconds. auto, the default, takes the slowest search of the run of
# the input (column 18 of the .perf next to it, in ms) minus its busiest thread; without them, the
# measured time on the first number of nodes minus the simulated one there; without that, 0.
# --warmup off keeps the first tasks as they are, as tasks like the others. --cores off gives all
# the cores the same speed, so that the time of each task is its size, as measured.
#
# --plot prefix writes the figures, with Plots.jl and StatsPlots.jl, without opening a window:
# prefix-speedup.pdf and .png, the measured against the hypothetical speedup; and prefix.pdf and
# .png, the same above the times of the tasks, the simulated time against the number of nodes, and
# its spread over the trials.

using Printf
using Random
using Statistics

function parse_args(args)
    options = Dict("nodes" => "", "layout" => "", "trials" => "10000", "seed" => "1",
                   "column" => "0", "plot" => "", "measured" => "", "results" => "",
                   "overhead" => "auto", "warmup" => "auto", "cores" => "auto")
    files = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if startswith(a, "--")
            key = a[3:end]
            haskey(options, key) || error("Unknown option $a")
            i < length(args) || error("Option $a needs a value")
            options[key] = args[i + 1]
            i += 2
        else
            push!(files, a)
            i += 1
        end
    end
    isempty(files) && error("Usage: simulate_tasks.jl times... [--layout PxT] [--nodes 1,2,4,8] " *
                            "[--trials N] [--seed S] [--column C] [--results table | " *
                            "--measured nodes:seconds,...] [--overhead auto|seconds] " *
                            "[--warmup auto|off] [--cores auto|off] [--plot prefix]")
    options["warmup"] in ("auto", "off") || error("--warmup is auto or off")
    options["cores"] in ("auto", "off") || error("--cores is auto or off")
    return (files = files,
            layout = isempty(options["layout"]) ? nothing : Tuple(parse.(Int, split(options["layout"], "x"))),
            nodes = parse.(Int, split(options["nodes"], ",", keepempty = false)),
            trials = parse(Int, options["trials"]),
            seed = parse(Int, options["seed"]),
            column = parse(Int, options["column"]),
            plot = options["plot"],
            results = options["results"],
            measured = [parse.(Float64, split(m, ":")) for m in split(options["measured"], ",", keepempty=false)],
            overhead = options["overhead"],
            warmup = options["warmup"] == "auto",
            cores = options["cores"] == "auto")
end

# The settings of one file of times: its layout, numbers of nodes, measurements and figures.
function settings(o, file, layout)
    processes, threads = o.layout !== nothing ? o.layout : layout
    measured = !isempty(o.results) ? read_results(o.results, "$(processes)x$(threads)") :
               Dict(Int(n) => t for (n, t) in o.measured)
    nodes = !isempty(o.nodes) ? o.nodes : !isempty(measured) ? sort(collect(keys(measured))) : [1, 2, 4, 8]
    plot = isempty(o.plot) || length(o.files) == 1 ? o.plot : "$(o.plot)-$(processes)x$(threads)"
    return (file = file, processes = processes, threads = threads, nodes = nodes,
            trials = o.trials, seed = o.seed, plot = plot, measured = measured,
            overhead = o.overhead, warmup = o.warmup, cores = o.cores)
end

# The measured time in seconds of each number of nodes, from the lines of layout in a table of
# 020-concat: the slowest search of a process (column 21, in ms), averaged over the trials.
function read_results(file, layout)
    times = Dict{Int, Vector{Float64}}()
    for line in eachline(file)
        fields = split(line, '\t')
        length(fields) >= 21 && fields[2] == layout || continue
        push!(get!(times, parse(Int, fields[1]), Float64[]), parse(Float64, fields[21]) / 1000)
    end
    isempty(times) && println("No measurement of the layout $layout in $file")
    return Dict(n => mean(t) for (n, t) in times)
end

# The times of the tasks, in the order of the code if the input gives it (ordered = true), else in
# that of the input; the thread of each, if the input gives them, else nothing; the layout of those
# threads, else that of the name of the file, else nothing; the busy time of the busiest thread, if
# known; and the slowest search of the run in seconds, from the .perf next to the file, if any.
function read_times(file, column)
    times, places, threads = Float64[], Int[], String[]
    for (n, line) in enumerate(eachline(file))
        line = strip(line)
        (isempty(line) || startswith(line, "#")) && continue
        fields = split(line)
        field = column == 0 ? fields[end] : fields[column]
        t = tryparse(Float64, field)
        if t === nothing
            isempty(times) && n == 1 && continue  # header
            error("$file:$n: not a time: $field")
        end
        t >= 0 || error("$file:$n: negative time: $t")
        push!(times, t)
        for f in fields
            (m = match(r"^conf(\d+)$", f)) === nothing || push!(places, parse(Int, m[1]))
            occursin(r"^p\d+t\d+$", f) && push!(threads, f)
        end
    end
    isempty(times) && error("$file: no times")

    has_threads = length(threads) == length(times)
    ordered = length(places) == length(times)
    if ordered
        sort(places) == 0:length(times)-1 || error("$file: the places conf0... are not those of all the tasks")
        times[places .+ 1] = copy(times)
        has_threads && (threads[places .+ 1] = copy(threads))
    end

    layout, busiest = nothing, nothing
    if has_threads
        pt = [parse.(Int, split(t[2:end], "t")) for t in threads]
        layout = (maximum(first.(pt)) + 1, maximum(last.(pt)) + 1)
        busy = Dict{String, Float64}()
        for (thread, t) in zip(threads, times)
            busy[thread] = get(busy, thread, 0.0) + t
        end
        busiest = maximum(values(busy))
    elseif (m = match(r"\.n\d+\.(\d+)x(\d+)\.", basename(file))) !== nothing
        layout = (parse(Int, m[1]), parse(Int, m[2]))
    end

    search = nothing
    perf = replace(file, r"\.tsv$" => ".perf")
    if perf != file && isfile(perf)
        fields = split(strip(read(perf, String)), '\t')
        length(fields) >= 18 && (search = parse(Float64, fields[18]) / 1000)
    end
    return (times = times, ordered = ordered, threads = has_threads ? threads : nothing,
            layout = layout, busiest = busiest, search = search)
end

# The sizes of the tasks, the speeds of the cores and the warm-ups of the threads (see the top of the
# file), from the times as measured:
# - times: the size of each task, its time divided by the speed of its thread, and for the first
#   task of each thread at most the median of the others;
# - speeds and pool: the speed and the warm-up of each thread, in the order of their first tasks
#   (the warm-up in the units of the sizes: a thread pays it times its speed);
# - slots and speed_slots: for the order of the code, the warm-up and speed of the thread that took
#   each of the first threads tasks of the block of each process, on 1 node with processes
#   processes (0 and 1 for a task that was not the first of its thread), as threads x processes
#   matrices;
# - first and factor: whether each task was the first of its thread, and the speed of its thread.
# Without the threads in the input, the sizes are the times, all the speeds 1, and no warm-up;
# warmup = false keeps the first tasks as they are; cores = false gives all the threads speed 1.
function split_times(input, processes, threads, warmup, cores)
    n = length(input.times)
    if input.threads === nothing
        return (times = input.times, pool = Float64[], speeds = Float64[],
                slots = zeros(threads, processes), speed_slots = ones(threads, processes),
                first = falses(n), factor = ones(n))
    end
    first_task = Dict{String, Int}()  # thread => its first task, in the order of the code
    for (c, thread) in enumerate(input.threads)
        haskey(first_task, thread) || (first_task[thread] = c)
    end
    is_first = falses(n)
    is_first[collect(values(first_task))] .= true

    # The speed of each thread: the median of its tasks after the first against that of all of them.
    later = Dict{String, Vector{Float64}}()
    for c in findall(.!is_first)
        push!(get!(later, input.threads[c], Float64[]), input.times[c])
    end
    typical = median(input.times[.!is_first])
    speed = Dict(thread => cores && length(get(later, thread, Float64[])) >= 3 ?
                           median(later[thread]) / typical : 1.0
                 for thread in keys(first_task))
    factor = [speed[thread] for thread in input.threads]

    times = input.times ./ factor
    warm = zeros(n)
    if warmup
        typical = median(times[.!is_first])
        for c in findall(is_first)
            warm[c] = max(0.0, times[c] - typical)
            times[c] -= warm[c]
        end
    end
    slots, speed_slots = zeros(threads, processes), ones(threads, processes)
    for p in 0:processes-1, i in 1:threads
        c = div(p * n, processes) + i
        if c <= div((p + 1) * n, processes)
            slots[i, p + 1], speed_slots[i, p + 1] = warm[c], factor[c]
        end
    end
    firsts = findall(is_first)
    return (times = times, pool = warmup ? warm[firsts] : Float64[],
            speeds = cores ? factor[firsts] : Float64[], slots = slots, speed_slots = speed_slots,
            first = is_first, factor = factor)
end

# Time of one process whose threads take the tasks of block in order: at the start one each, with
# its warm-up, then each the next one as soon as it is free. A task takes its size times the speed
# of its thread, in speed. free holds the time at which each thread is free.
function process_time!(free, block, warm, speed)
    fill!(free, 0.0)
    for (i, t) in enumerate(block)
        if i <= length(free)
            free[i] = (warm[i] + t) * speed[i]
        else
            j = argmin(free)
            free[j] += t * speed[j]
        end
    end
    return maximum(free)
end

# Time of the tasks in order: the slowest of the num_parts processes, the warm-ups and speeds of the
# threads of process part in column part + 1 of warm and speed.
function order_time!(order, free, times, num_parts, warm, speed)
    n = length(order)
    slowest = 0.0
    for part in 0:num_parts-1
        # The block of part, as in shuffled_configurations.
        first, last = div(part * n, num_parts) + 1, div((part + 1) * n, num_parts)
        slowest = max(slowest, process_time!(free, @view(times[order[first:last]]),
                                             @view(warm[:, part + 1]), @view(speed[:, part + 1])))
    end
    return slowest
end

# Draws the warm-ups and speeds of the threads of num_parts processes, with processes processes per
# node, into warm and speed: on each node, those measured, of the same threads, in a random order, if
# there are as many as threads on a node, else drawn with replacement; no warm-up and speed 1 if none
# were measured.
function draw_threads!(warm, speed, pool, speeds, processes, rng)
    threads, num_parts = size(warm)
    n = max(length(pool), length(speeds))
    fill!(warm, 0.0)
    fill!(speed, 1.0)
    n == 0 && return
    per_node = threads * processes
    for node in 0:div(num_parts, processes)-1
        pick = n == per_node ? randperm(rng, n) : rand(rng, 1:n, per_node)
        columns = node * processes + 1:(node + 1) * processes
        isempty(pool) || (warm[:, columns] = reshape(pool[pick], threads, processes))
        isempty(speeds) || (speed[:, columns] = reshape(speeds[pick], threads, processes))
    end
end

# The times of trials random runs on num_parts processes, processes per node. Trial i shuffles the
# tasks in the same order for any num_parts, with its own random generator, and draws the warm-ups
# and speeds with another.
function simulate(times, pool, speeds, num_parts, processes, threads, trials, seed)
    results = Vector{Float64}(undef, trials)
    chunks = collect(Iterators.partition(1:trials, cld(trials, Threads.nthreads())))
    Threads.@threads for c in eachindex(chunks)
        order_rng, thread_rng = Xoshiro(seed + c), Xoshiro(seed + c + 1_000_000)
        order = collect(eachindex(times))
        free = zeros(threads)
        warm, speed = zeros(threads, num_parts), ones(threads, num_parts)
        for i in chunks[c]
            shuffle!(order_rng, order)
            draw_threads!(warm, speed, pool, speeds, processes, thread_rng)
            results[i] = order_time!(order, free, times, num_parts, warm, speed)
        end
    end
    return results
end

function main(args)
    o = parse_args(args)
    for (i, file) in enumerate(o.files)
        i > 1 && println("\n", "-"^110, "\n")
        input = read_times(file, o.column)
        o.layout === nothing && input.layout === nothing &&
            error("$file: no layout; give one with --layout")
        run(settings(o, file, input.layout), input)
    end
end

function run(o, input)
    threads = o.processes * o.threads  # per node
    split = split_times(input, o.processes, o.threads, o.warmup, o.cores)
    times = split.times
    total, longest = sum(times), maximum(times)
    warmup = !isempty(split.pool)
    cores = !isempty(split.speeds)

    println(o.file)
    @printf("%d tasks as measured: sum %.4g s, mean %.4g, median %.4g, max %.4g, coefficient of variation %.3f\n",
            length(times), sum(input.times), mean(input.times), median(input.times),
            maximum(input.times), std(input.times) / mean(input.times))
    if cores
        slow = sort([(s, t) for (s, t) in zip(split.speeds, input.threads[split.first]) if s >= 1.05], rev = true)
        @printf("speeds of the %d threads, against the median: median %.3f, IQR %.3f–%.3f; %d at least 5%% slower: %s\n",
                length(split.speeds), median(split.speeds), quantile(split.speeds, 0.25),
                quantile(split.speeds, 0.75), length(slow),
                isempty(slow) ? "none" : join([@sprintf("%s %.2f", t, s) for (s, t) in slow], ", "))
    end
    @printf("%d tasks as sizes%s: mean %.4g s, median %.4g, max %.4g, coefficient of variation %.3f\n",
            length(times), cores && warmup ? ", without the speeds of the cores and the warm-up" :
                           cores ? ", without the speeds of the cores" : warmup ? ", without the warm-up" : "",
            mean(times), median(times), longest, std(times) / mean(times))
    if warmup
        @printf("warm-up of the %d threads, beyond the median task: mean %.4g s, median %.4g, max %.4g\n",
                length(split.pool), mean(split.pool), median(split.pool), maximum(split.pool))
    end
    @printf("layout %dx%d per node, %d trials per number of nodes, %d Julia threads\n",
            o.processes, o.threads, o.trials, Threads.nthreads())

    # The tasks alone, then on the cores as measured, then with the warm-ups, all before the overhead.
    sim(pool, speeds) = [simulate(times, pool, speeds, n * o.processes, o.processes, o.threads, o.trials, o.seed)
                         for n in o.nodes]
    alone = sim(Float64[], Float64[])
    cored = cores ? sim(Float64[], split.speeds) : alone
    warmed = warmup ? sim(split.pool, split.speeds) : cored
    # All the tasks of equal length, on cores of equal speed: no order and no trials matter.
    equal = [order_time!(collect(eachindex(times)), zeros(o.threads), fill(mean(times), length(times)),
                         n * o.processes, zeros(o.threads, n * o.processes), ones(o.threads, n * o.processes))
             for n in o.nodes]
    fixed = input.ordered ?
            [order_time!(collect(eachindex(times)), zeros(o.threads), times, n * o.processes,
                         repeat(split.slots, 1, n), repeat(split.speed_slots, 1, n)) for n in o.nodes] : nothing
    if fixed !== nothing && input.busiest !== nothing && o.nodes[1] == 1
        @printf("check: the order of the code on 1 node takes %.5g s simulated, and its busiest thread %.5g s measured\n",
                fixed[1], input.busiest)
        # Two threads that finish within the gap between tasks of the run may take the next task in
        # the other order than in the run; with threads of different speeds, the tasks after it then
        # go to threads of other speeds than in the run.
        cores && abs(fixed[1] - input.busiest) > 1e-3 * input.busiest &&
            println("       (the speeds of the cores differ, so threads that finished almost together and took " *
                    "the next tasks in the other order than in the run give other times)")
    end

    overhead = o.overhead == "auto" ? nothing : parse(Float64, o.overhead)
    if overhead !== nothing
        @printf("overhead: %.4g s\n", overhead)
    elseif input.search !== nothing && input.busiest !== nothing
        overhead = input.search - input.busiest
        @printf("overhead: %.4g s, the slowest search of the run, %.5g s, minus its busiest thread\n",
                overhead, input.search)
    elseif haskey(o.measured, o.nodes[1])
        overhead = o.measured[o.nodes[1]] - (fixed === nothing ? mean(warmed[1]) : fixed[1])
        @printf("overhead: %.4g s, the measurement on %d node(s) minus the simulation there\n",
                overhead, o.nodes[1])
    else
        overhead = 0.0
    end
    overhead < 0 && println("The overhead is negative: the measurements disagree with the simulation")
    println()

    @printf("%5s %9s %9s %10s %9s %10s %10s %10s %10s %9s %9s\n", "nodes", "processes", "threads",
            "mean", "95% ci +-", "median", "p95", "code order", "ideal", "speedup", "ideal")
    rows = []
    for (k, nodes) in enumerate(o.nodes)
        num_parts = nodes * o.processes
        r = warmed[k] .+ overhead
        m = mean(r)
        ci = 1.96 * std(r) / sqrt(length(r))
        # No schedule does better than an even share of the work over the speeds of the cores, as
        # every thread takes a task, or than the longest task on the fastest core.
        capacity = cores ? threads * mean(1 ./ split.speeds) : threads
        share = total / (nodes * capacity) + (warmup ? mean(split.pool) : 0.0)
        ideal = max(share, longest * (cores ? minimum(split.speeds) : 1.0)) + overhead
        f = fixed === nothing ? NaN : fixed[k] + overhead
        push!(rows, (nodes = nodes, results = r, mean = m, ideal = ideal, fixed = f,
                     equal = equal[k], alone = mean(alone[k]), cored = mean(cored[k]),
                     warmed = mean(warmed[k])))
        @printf("%5d %9d %9d %10.5g %9.3g %10.5g %10.5g %10.5g %10.5g %9.2f %9.2f\n",
                nodes, num_parts, nodes * threads, m, ci, median(r), quantile(r, 0.95), f,
                ideal, rows[1].mean / m, rows[1].ideal / ideal)
    end
    println("\ntimes in s, of the tasks, speeds of the cores, warm-up and overhead; code order: in the order of the code; " *
            "speedup: of the mean,\nagainst the first row; ideal: of the lower bound " *
            "max(even share, longest task) + overhead")

    comparison = compare_speedups(o, rows)
    if !isempty(comparison)
        @printf("\n%5s %9s %8s %8s %7s %7s %8s %9s %10s %12s %11s %10s %10s\n", "nodes", "measured",
                "speedup", "perfect", "equal", "tasks", "+ cores", "+ warm-up", "+ overhead", "5–95%",
                "code order", "explained", "code order")
        percent(x) = isnan(x) ? "–" : @sprintf("%.0f%%", 100 * x)
        for c in comparison
            @printf("%5d %9.5g %8.2f %8.2f %7.2f %7.2f %8.2f %9.2f %10.2f %5.2f–%-6.2f %11.2f %10s %10s\n",
                    c.nodes, c.time, c.measured, c.perfect, c.equal, c.alone, c.cored, c.warmed,
                    c.hypothetical, c.low, c.high, c.fixed, percent(c.explained),
                    percent(c.explained_fixed))
        end
        println("\nspeedups against $(rows[1].nodes) node(s), of the mean over random orders unless " *
                "said; measured: in s, and its speedup;\nequal: of tasks of equal length; tasks: of the sizes of the tasks alone; " *
                "+ cores: on cores of the speeds measured; + warm-up: and the\nwarm-ups; + overhead: and the overhead, the hypothetical speedup, with the range of " *
                "single runs; code order: the same, in the order of the code; explained: (perfect - " *
                "hypothetical) /\n(perfect - measured), of the mean and of the order of the code")
    end

    if !isempty(o.plot)
        ENV["GKSwstype"] = "nul"  # GR draws into the files only, and opens no window
        @eval using Plots, StatsPlots
        Base.invokelatest(plot_results, o, input, split, rows, comparison, overhead)
    end
end

# The measured speedups next to the hypothetical ones, for the numbers of nodes that were both
# measured and simulated; none without a measurement on the first number of nodes. The range of
# the hypothetical speedup is that of single runs (see hypothetical_speedups).
function compare_speedups(o, rows)
    measured = o.measured
    base = rows[1]
    haskey(measured, base.nodes) || return []
    for n in keys(measured)
        any(r -> r.nodes == n, rows) || println("No simulation for the measurement on $n nodes")
    end
    comparison = []
    for r in rows
        haskey(measured, r.nodes) || continue
        s = measured[base.nodes] / measured[r.nodes]
        h = base.mean / r.mean
        ratios = base.results ./ r.results  # as in hypothetical_speedups
        perfect = r.nodes / base.nodes
        f = base.fixed / r.fixed
        lost = perfect - s
        explained(x) = lost > 1e-9 && !isnan(x) ? (perfect - x) / lost : NaN
        push!(comparison, (nodes = r.nodes, time = measured[r.nodes], measured = s,
                           equal = base.equal / r.equal, alone = base.alone / r.alone,
                           cored = base.cored / r.cored, warmed = base.warmed / r.warmed,
                           hypothetical = h, low = quantile(ratios, 0.05),
                           high = quantile(ratios, 0.95), fixed = f, perfect = perfect,
                           explained = explained(h), explained_fixed = explained(f)))
    end
    return comparison
end

# The colors of the figures: the simulation, the measurements, the simulation of the order of the
# code, and the references in gray.
const SIMULATED = "#2a78d6"
const MEASURED = "#eb6834"
const FIXED = "#1baf7a"
const WARMUP = "#6250d6"
const SLOW = "#e87ba4"
const REFERENCE = "#52514e"
const MUTED = "#898781"
const GRID = "#e1e0d9"

# Ticks at 1, 2 and 5 times the powers of 10 in the range of values, with plain labels; on a range
# of less than a factor of 5, at more steps of each power of 10.
function nice_ticks(values)
    lo, hi = extrema(filter(>(0), values))
    steps = hi / lo < 5 ? (1:9) : hi / lo < 50 ? (1, 2, 3, 5) : (1, 2, 5)
    ticks = [s * 10.0^e for e in floor(Int, log10(lo)):ceil(Int, log10(hi)) for s in steps]
    ticks = filter(t -> lo / 1.05 <= t <= hi * 1.05, ticks)
    return (ticks, [t >= 1 ? @sprintf("%.0f", t) : @sprintf("%.3g", t) for t in ticks])
end

function plot_results(o, input, split, rows, comparison, overhead)
    times = split.times
    gr()
    default(fontfamily = "Helvetica", titlefontsize = 11, guidefontsize = 9, tickfontsize = 8,
            legendfontsize = 8, linewidth = 2, framestyle = :axes, grid = :y, gridcolor = GRID,
            gridalpha = 1, gridlinewidth = 0.5, foreground_color_axis = MUTED,
            foreground_color_border = MUTED, foreground_color_legend = nothing,
            background_color_legend = nothing, titlelocation = :left)

    nodes = [r.nodes for r in rows]
    means = [r.mean for r in rows]
    low = [quantile(r.results, 0.05) for r in rows]
    high = [quantile(r.results, 0.95) for r in rows]
    bound = [r.ideal for r in rows]
    perfect = means[1] .* nodes[1] ./ nodes
    total, longest = sum(times), maximum(times)
    layout = "$(o.processes)x$(o.threads) per node"
    plus = overhead > 0 ? @sprintf(" + %.0f s overhead", overhead) : ""
    measured_nodes = sort(collect(keys(o.measured)))
    measured_times = [o.measured[n] for n in measured_nodes]
    node_ticks = (nodes, string.(nodes))

    # The tasks as measured: in the order of the code, with the first task of each thread and the
    # tasks on slow cores apart, if known; else longest first.
    if input.ordered && input.threads !== nothing
        firsts = findall(split.first)
        slow = findall(.!split.first .& (split.factor .>= 1.05))
        rest = findall(.!split.first .& (split.factor .< 1.05))
        tasks = scatter(rest, input.times[rest], color = SIMULATED, markersize = 2.5,
                        markerstrokewidth = 0, label = "task",
                        xlabel = "task, in the order of the code", ylabel = "time (s)",
                        legend = :topright, title = "")
        isempty(slow) || scatter!(tasks, slow, input.times[slow], color = SLOW, markersize = 3,
                                  markerstrokewidth = 0, label = "task on a core at least 5% slower")
        scatter!(tasks, firsts, input.times[firsts], color = WARMUP, markersize = 3.5,
                 markerstrokewidth = 0, label = "first task of a thread, with its warm-up")
        typical = median(input.times[.!split.first])
        hline!(tasks, [typical], color = REFERENCE, linewidth = 1, linestyle = :dash,
               label = @sprintf("median of the others, %.0f s", typical))
        ylims!(tasks, (0, maximum(input.times) * 1.25))
    else
        sorted = sort(times, rev = true)
        log_y = minimum(sorted) > 0 && maximum(sorted) / minimum(sorted) > 20
        tasks = plot(1:length(sorted), sorted, color = SIMULATED, label = "task",
                     yscale = log_y ? :log10 : :identity,
                     yticks = log_y ? nice_ticks(sorted) : :auto, legend = :bottomleft,
                     xlabel = "task, longest first", ylabel = "time (s)", title = "")
        for (r, style) in ((rows[1], :dash), (rows[end], :dot))
            share = total / (r.nodes * o.processes * o.threads)
            hline!(tasks, [share], color = REFERENCE, linewidth = 1, linestyle = style,
                   label = "even share per thread, $(r.nodes) node$(r.nodes == 1 ? "" : "s")")
        end
    end

    # The time against the number of nodes, on log scales, where perfect scaling is a line.
    elapsed = plot(nodes, means, ribbon = (means .- low, high .- means), fillalpha = 0.15,
                   color = SIMULATED, marker = :circle, markersize = 4, markerstrokewidth = 0,
                   label = "simulated, mean and 5–95%", xscale = :log2, yscale = :log10,
                   xticks = node_ticks, legend = :bottomleft,
                   yticks = nice_ticks(vcat(low, high, perfect, bound, measured_times)),
                   xlabel = "nodes", ylabel = "time (s)", title = "")
    plot!(elapsed, nodes, perfect, color = REFERENCE, linewidth = 1, linestyle = :dash,
          label = "perfect scaling")
    plot!(elapsed, nodes, bound, color = REFERENCE, linewidth = 1, linestyle = :dot,
          label = "lower bound")
    hline!(elapsed, [longest + overhead], color = MUTED, linewidth = 0.5, label = "longest task" * plus)
    if input.ordered
        scatter!(elapsed, nodes, [r.fixed for r in rows], color = FIXED, marker = :diamond,
                 markersize = 5, markerstrokewidth = 0, label = "simulated, order of the code")
    end
    if !isempty(o.measured)
        scatter!(elapsed, measured_nodes, measured_times, color = MEASURED, markersize = 5,
                 markerstrokewidth = 0, label = "measured")
    end
    headline = speedup_plots(o, input, split, rows, comparison, overhead, "Measured and hypothetical speedup, $layout")

    # The spread of the simulated time over the trials, for each number of nodes.
    spread = boxplot(repeat(string.(nodes), inner = o.trials), reduce(vcat, [r.results for r in rows]),
                     color = SIMULATED, fillalpha = 0.25, linecolor = SIMULATED, linewidth = 1,
                     outliers = false, whisker_range = Inf, label = "",
                     xlabel = "nodes", ylabel = "time (s)", title = "")
    hline!(spread, [longest + overhead], color = MUTED, linewidth = 0.5, label = "longest task" * plus)

    title!(tasks, "b. Times of the $(length(times)) tasks, on 1 node")
    title!(elapsed, "c. Time of the search")
    title!(spread, "d. Spread over $(o.trials) trials")
    # The speedup above, the other three in a row below it.
    grid_layout = grid(2, 1, heights = [0.55, 0.45])
    grid_layout[1, 1] = grid(1, 2, widths = [0.55, 0.45])
    grid_layout[2, 1] = grid(1, 3)
    overview = plot(speedup_plots(o, input, split, rows, comparison, overhead, "a. Speedup")..., tasks, elapsed, spread,
                    layout = grid_layout, size = (1300, 1000),
                    left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, top_margin = 2Plots.mm,
                    plot_title = "$(replace(basename(o.file), r"\.n\d+\..*" => "")): " *
                                 "$(length(times)) tasks, $(layout)$(plus)",
                    plot_titlefontsize = 12, plot_titlelocation = :left)
    alone = plot(headline..., layout = grid(1, 2, widths = [0.55, 0.45]), size = (1300, 600),
                 left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, top_margin = 2Plots.mm)
    for (figure, name) in ((alone, "$(o.plot)-speedup"), (overview, o.plot)),
        extension in ("pdf", "png")
        savefig(figure, "$(name).$(extension)")
        println("wrote $(name).$(extension)")
    end
end

# The measured speedup against the hypothetical one: next to perfect scaling, on linear scales,
# where it is the diagonal and the gap between a line and the diagonal is the speedup lost; and
# close up, on the range of the two, each labelled with its value at the last number of nodes.
function speedup_plots(o, input, split, rows, comparison, overhead, title)
    return (speedup_plot(o, input, split, rows, comparison, overhead, title),
            speedup_closeup(o, input, rows, comparison, "Close up"))
end

# The hypothetical speedup of each number of nodes, and its 5–95% range over single runs: trial i on
# the first number of nodes against trial i on the others, which shuffle the tasks in the same order,
# as the code does on any number of nodes.
function hypothetical_speedups(rows)
    base = rows[1]
    ratios = [base.results ./ r.results for r in rows]
    return (base.mean ./ [r.mean for r in rows], quantile.(ratios, 0.05), quantile.(ratios, 0.95))
end

function speedup_closeup(o, input, rows, comparison, title)
    nodes = [r.nodes for r in rows]
    hypothetical, low, high = hypothetical_speedups(rows)
    fixed = rows[1].fixed ./ [r.fixed for r in rows]
    measured = [c.measured for c in comparison]
    lo, hi = extrema(vcat(low, high, measured, input.ordered ? fixed : Float64[]))
    pad = 0.08 * max(hi - lo, 0.1)
    p = plot(nodes, hypothetical, ribbon = (hypothetical .- low, high .- hypothetical),
             fillalpha = 0.15, color = SIMULATED, linewidth = 2.5, marker = :circle, markersize = 6,
             markerstrokewidth = 0, label = "", xticks = (nodes, string.(nodes)),
             xlims = (0, maximum(nodes) * 1.18), ylims = (lo - pad, hi + pad),
             xlabel = "nodes", ylabel = "speedup", guidefontsize = 10, tickfontsize = 9,
             titlefontsize = 13, title = title)
    # Perfect scaling, where it is in range: only its first points, at small speedups.
    plot!(p, nodes, nodes ./ nodes[1], color = REFERENCE, linewidth = 1, linestyle = :dash, label = "")
    labels = [(hypothetical[end], SIMULATED, "hypothetical")]
    if input.ordered
        scatter!(p, nodes, fixed, color = FIXED, marker = :diamond, markersize = 7,
                 markerstrokewidth = 0, label = "")
    end
    if !isempty(comparison)
        plot!(p, [c.nodes for c in comparison], measured, color = MEASURED, linewidth = 2.5,
              marker = :circle, markersize = 6, markerstrokewidth = 0, label = "")
        comparison[end].nodes == nodes[end] && push!(labels, (measured[end], MEASURED, "measured"))
    end
    # The higher of the two labels above its point, the lower below it, so that they do not meet.
    top = argmax(first.(labels))
    for (i, (value, color, _)) in enumerate(labels)
        annotate!(p, nodes[end] * 1.025, value,
                  text(@sprintf("%.2f×", value), 10, color, :left,
                       length(labels) == 1 ? :vcenter : i == top ? :bottom : :top))
    end
    return p
end

function speedup_plot(o, input, split, rows, comparison, overhead, title)
    nodes = [r.nodes for r in rows]
    base = rows[1]
    unit = "$(base.nodes) node$(base.nodes == 1 ? "" : "s")"
    hypothetical, low, high = hypothetical_speedups(rows)
    perfect = nodes ./ base.nodes
    equal = base.equal ./ [r.equal for r in rows]
    alone = base.alone ./ [r.alone for r in rows]
    cored = base.cored ./ [r.cored for r in rows]
    fixed = base.fixed ./ [r.fixed for r in rows]
    # Each entry of the legend gives its speedup on the last number of nodes.
    at(x) = @sprintf("%.2f×", x[end])

    p = plot(nodes, perfect, color = REFERENCE, linewidth = 1, linestyle = :dash,
             label = "perfect scaling: $(at(perfect)) on $(nodes[end]) nodes",
             xticks = (nodes, string.(nodes)),
             xlims = (0, maximum(nodes) * 1.04), ylims = (0, maximum(perfect) * 1.04),
             xlabel = "nodes", ylabel = "speedup over $unit", legend = :topleft,
             legendfontsize = 9, guidefontsize = 10, tickfontsize = 9, titlefontsize = 13,
             title = title)
    plot!(p, nodes, equal, color = REFERENCE, linewidth = 1, linestyle = :dashdot,
          label = "tasks of equal length: $(at(equal))")
    cores = !isempty(split.speeds)
    if cores || !isempty(split.pool) || overhead > 0
        plot!(p, nodes, alone, color = SIMULATED, linewidth = 1.5, linestyle = :dot,
              label = "tasks alone" * (cores ? ", on cores of equal speed" : "") * ": $(at(alone))")
    end
    if cores && (!isempty(split.pool) || overhead > 0)
        plot!(p, nodes, cored, color = SLOW, linewidth = 1.5, linestyle = :dot,
              label = "tasks on cores of the speeds measured: $(at(cored))")
    end
    plot!(p, nodes, hypothetical, ribbon = (hypothetical .- low, high .- hypothetical),
          fillalpha = 0.15, color = SIMULATED, linewidth = 2.5, marker = :circle, markersize = 6,
          markerstrokewidth = 0,
          label = "hypothetical: tasks" * (cores ? " + cores" : "") *
                  (isempty(split.pool) ? "" : " + warm-up") *
                  (overhead > 0 ? " + overhead" : "") * ": $(at(hypothetical)) (5–95% of runs)")
    if input.ordered
        scatter!(p, nodes, fixed, color = FIXED, marker = :diamond, markersize = 7,
                 markerstrokewidth = 0, label = "the same, in the order of the code: $(at(fixed))")
    end
    if !isempty(comparison)
        measured = [c.measured for c in comparison]
        plot!(p, [c.nodes for c in comparison], measured, color = MEASURED, linewidth = 2.5,
              marker = :circle, markersize = 6, markerstrokewidth = 0,
              label = "measured: $(at(measured))")
        # One short line under the legend, left of the diagonal.
        last = comparison[end]
        if !isnan(last.explained)
            annotate!(p, 0.03 * maximum(nodes), 0.7 * maximum(perfect),
                      text(@sprintf("The hypothetical explains %.0f%% of the speedup lost", 100 * last.explained) *
                           (isnan(last.explained_fixed) ? "" :
                            @sprintf("\n(in the order of the code, %.0f%%)", 100 * last.explained_fixed)),
                           10, REFERENCE, :left, :top))
        end
    elseif !isempty(o.measured)
        annotate!(p, 0.03 * maximum(nodes), 0.7 * maximum(perfect),
                  text("No measurement on $unit to compare with", 10, REFERENCE, :left, :top))
    end
    return p
end

main(ARGS)
