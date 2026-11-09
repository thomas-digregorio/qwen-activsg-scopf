# Download the TAMU ACTIVSg25k dataset from Google Drive.
# Drive file ID from https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg25k/
$ErrorActionPreference = "Stop"
$id  = "1r5RIJefM3q2Bb59-FEA62h3yfMC-3ORg"
$out = "C:\Users\thoma\Documents\OpenCodeProjects\OPF-scratch\activsg25k_dataset.zip"
$tmp = "C:\Users\thoma\AppData\Local\Temp\opencode\activsg25k_probe.html"

# First request: may be the file itself, or a virus-scan / confirm page
curl.exe -sL -o "$tmp" "https://drive.google.com/uc?export=download&id=$id"
$head = Get-Content "$tmp" -TotalCount 1 -ErrorAction SilentlyContinue
Write-Output "first-response head: $($head.Substring(0,[Math]::Min(60,$head.Length)))"

# If it's an HTML confirm page, extract the token and retry
if ($head -like '<!DOCTYPE*' -or $head -like '<html*' -or $head -like '<*') {
    $html = Get-Content "$tmp" -Raw
    $token = $null
    $m = [regex]::Match($html, 'name="confirm" value="([^"]+)"')
    if ($m.Success) { $token = $m.Groups[1].Value }
    if (-not $token) {
        $m2 = [regex]::Match($html, 'confirm=([0-9A-Za-z_-]+)')
        if ($m2.Success) { $token = $m2.Groups[1].Value }
    }
    if (-not $token) {
        $m3 = [regex]::Match($html, 'href="(https://drive\.google\.com/uc[^"]*export=download[^"]*)"')
        if ($m3.Success) {
            curl.exe -sL -o "$out" $m3.Groups[1].Value
            Write-Output "downloaded via embedded link"
        } else {
            throw "NO DOWNLOAD LINK FOUND; page saved to $tmp"
        }
    } else {
        curl.exe -sL -o "$out" "https://drive.google.com/uc?export=download&confirm=$token&id=$id"
        Write-Output "downloaded with confirm token"
    }
} else {
    Move-Item -Force "$tmp" "$out"
    Write-Output "downloaded directly"
}
Get-Item $out | Select-Object Name, @{n='MB';e={[math]::Round($_.Length/1MB,2)}}
# quick integrity: must be a ZIP (PK signature)
$b = [System.IO.File]::ReadAllBytes($out)[0..1]
if ($b[0] -eq 0x50 -and $b[1] -eq 0x4B) { "ZIP signature OK" } else { "NOT A ZIP - inspect $out" }
