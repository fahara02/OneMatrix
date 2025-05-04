# PowerShell build script for the OneMatrix library
param(
    [switch]$clean,
    [switch]$run
)

# Configuration
$buildDir = "build"
$compiler = "clang++"
$src = "src/main.cpp"
$output = "$buildDir/matrix_test.exe"
$flags = "-Wall -Wextra -std=c++17 -fcolor-diagnostics -I./src"

# Create build directory if it doesn't exist
if (-not (Test-Path $buildDir)) {
    New-Item -ItemType Directory -Path $buildDir | Out-Null
    Write-Host "Created build directory: $buildDir"
}

# Clean build directory if requested
if ($clean) {
    Get-ChildItem -Path $buildDir -Recurse | Remove-Item -Force -Recurse
    Write-Host "Cleaned build directory"
    if (-not $run) { exit 0 }
}

# Compile
Write-Host "Compiling $src..."
$compileCommand = "$compiler $flags $src -o $output"
Write-Host "> $compileCommand"
Invoke-Expression $compileCommand

# Check if compilation was successful
if ($LASTEXITCODE -eq 0) {
    Write-Host "Compilation successful: $output"
    
    # Run the executable if requested
    if ($run) {
        Write-Host "Running $output..."
        Write-Host "Press Enter when done to close this window..." -ForegroundColor Yellow
        & $output
        Read-Host "`nExecution complete! Press Enter to exit"
    }
} else {
    Write-Host "Compilation failed with exit code $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}
