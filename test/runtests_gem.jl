using Test
using ProgressMeter

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Directories inside 'test/' to scan
const TEST_FOLDERS = ["unit", "examples"]

# The regex pattern provided
const FILE_PATTERN = r"^\d{2}(\.?\d?)*?_(?<testname>[A-Z][^_]*?)\.jl"

# Files to exclude (filenames only)
const EXCLUDE = [
# "03_BrokenTest.jl",
]

# If not empty, ONLY these files will run (filenames only), ignoring EXCLUDE
const ONLY = [
# "01_SpecificTest.jl",
]

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

struct TestFile
    path::String
    filename::String
    name::String
end

function get_test_files(folders, pattern)
    candidates = TestFile[]

    for folder in folders
        # Construct absolute path relative to this script
        folder_path = joinpath(@__DIR__, folder)

        # Skip if folder doesn't exist to prevent errors
        if !isdir(folder_path)
            @warn "Test folder not found: $folder_path"
            continue
        end

        for file in readdir(folder_path)
            m = match(pattern, file)
            if m !== nothing
                push!(candidates, TestFile(
                    joinpath(folder_path, file), # Full path for include()
                    file,                        # Filename for filtering
                    m[:testname]                 # Captured group for @testset name
                ))
            end
        end
    end
    return candidates
end

function filter_tests(candidates, only_list, exclude_list)
    if !isempty(only_list)
        return filter(t -> t.filename in only_list, candidates)
    else
        return filter(t -> !(t.filename in exclude_list), candidates)
    end
end

# ==============================================================================
# MAIN RUNNER
# ==============================================================================

@testset "Package Test Suite" begin
    # 1. Scan folders
    all_candidates = get_test_files(TEST_FOLDERS, FILE_PATTERN)

    # 2. Apply filtering logic
    tests_to_run = filter_tests(all_candidates, ONLY, EXCLUDE)

    if isempty(tests_to_run)
        @warn "No tests matched configuration!"
    else
        @info "Running $(length(tests_to_run)) test files..."

        # 3. Initialize Progress Bar
        p = Progress(length(tests_to_run); desc="Running Tests: ", color=:green)

        # 4. Execute
        for test_file in tests_to_run
            # Update progress bar with the name of the current test being run
            next!(p; showvalues=[(:Current, test_file.name)])

            @testset "$(test_file.name)" begin
                # Run the script in its own module scope if preferred, 
                # but standard practice is usually direct include
                include(test_file.path)
            end
        end
    end
end