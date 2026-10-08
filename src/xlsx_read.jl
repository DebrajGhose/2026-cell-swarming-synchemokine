using CodecZlib

function zip_entries(data::Vector{UInt8})
    u16(o) = UInt16(data[o+1]) | (UInt16(data[o+2]) << 8)
    u32(o) = UInt32(data[o+1]) | (UInt32(data[o+2]) << 8) |
             (UInt32(data[o+3]) << 16) | (UInt32(data[o+4]) << 24)

    eocd = -1
    for ii in (length(data) - 22):-1:max(0, length(data) - 66000)
        if u32(ii) == 0x06054b50; eocd = ii; break; end
    end
    eocd < 0 && error("not a zip file: no end-of-central-directory record")
    n  = Int(u16(eocd + 10))
    cd = Int(u32(eocd + 16))
    out = Dict{String,NamedTuple{(:method, :csize, :offset), Tuple{Int,Int,Int}}}()
    o = cd
    for ii in 1:n
        u32(o) == 0x02014b50 || error("bad central directory entry")
        method  = Int(u16(o + 10))
        csize   = Int(u32(o + 20))
        namelen = Int(u16(o + 28))
        extralen = Int(u16(o + 30))
        cmtlen  = Int(u16(o + 32))
        lo      = Int(u32(o + 42))
        name    = String(data[(o + 46 + 1):(o + 46 + namelen)])
        out[name] = (method = method, csize = csize, offset = lo)
        o += 46 + namelen + extralen + cmtlen
    end
    out
end

function zip_read(data::Vector{UInt8}, entries, name::AbstractString)
    haskey(entries, name) || error("no entry \"$name\" in the workbook")
    e = entries[name]
    u16(o) = UInt16(data[o+1]) | (UInt16(data[o+2]) << 8)
    o = e.offset
    UInt32(data[o+1]) | (UInt32(data[o+2]) << 8) | (UInt32(data[o+3]) << 16) |
        (UInt32(data[o+4]) << 24) == 0x04034b50 || error("bad local header for $name")
    dstart = o + 30 + Int(u16(o + 26)) + Int(u16(o + 28))
    raw = data[(dstart + 1):(dstart + e.csize)]
    e.method == 0 && return String(raw)
    e.method == 8 || error("unsupported zip compression method $(e.method) for $name")
    String(transcode(DeflateDecompressor, raw))
end

colof(ref) = (at = findfirst(isdigit, ref); at === nothing ? ref : ref[1:prevind(ref, at)])

function read_xlsx(path::AbstractString; want::Vector{String})
    data = read(path)
    ent  = zip_entries(data)

    shared = String[]
    if haskey(ent, "xl/sharedStrings.xml")
        s = zip_read(data, ent, "xl/sharedStrings.xml")
        for ii in eachmatch(r"<si>(.*?)</si>"s, s)
            push!(shared,
                  join(iii.captures[1] for iii in eachmatch(r"<t[^>]*>(.*?)</t>"s,
                                                            ii.captures[1])))
        end
    end

    wb   = zip_read(data, ent, "xl/workbook.xml")
    rels = zip_read(data, ent, "xl/_rels/workbook.xml.rels")
    target = Dict(ii.captures[1] => ii.captures[2] for ii in
                  eachmatch(r"Id=\"(rId\d+)\"[^>]*Target=\"(worksheets/sheet\d+\.xml)\"", rels))
    sheets = [(ii.captures[1], target[ii.captures[2]]) for ii in
              eachmatch(r"<sheet name=\"([^\"]*)\"[^>]*r:id=\"(rId\d+)\"", wb)]

    out = NamedTuple{(:name, :cols), Tuple{String, Dict{String,Vector{Float64}}}}[]
    for (sname, part) in sheets
        xml = zip_read(data, ent, "xl/" * part)
        rows = collect(eachmatch(r"<row r=\"(\d+)\"[^>]*>(.*?)</row>"s, xml))
        isempty(rows) && continue

        head = Dict{String,String}()
        for iii in eachmatch(r"<c r=\"([A-Z]+\d+)\"([^>]*)>(?:<v>([^<]*)</v>)?",
                             rows[1].captures[2])
            v = iii.captures[3]
            v === nothing && continue
            txt = occursin("t=\"s\"", iii.captures[2]) ? shared[parse(Int, v) + 1] : v
            head[colof(iii.captures[1])] = txt
        end
        wanted = Dict{String,String}()
        for (cl, iii) in head
            iii in want && (wanted[cl] = iii)
        end

        vals = Dict(iii => Float64[] for iii in want)
        for rowmatch in rows[2:end]
            got = Dict{String,Float64}()
            for iiii in eachmatch(r"<c r=\"([A-Z]+)\d+\"[^>]*>(?:<v>([^<]*)</v>)?",
                                  rowmatch.captures[2])
                cl = iiii.captures[1]
                haskey(wanted, cl) || continue
                v = iiii.captures[2]
                got[wanted[cl]] = (v === nothing || isempty(v)) ? NaN :
                                  something(tryparse(Float64, v), NaN)
            end
            for iiii in want
                push!(vals[iiii], get(got, iiii, NaN))
            end
        end
        push!(out, (name = sname, cols = vals))
    end
    out
end
