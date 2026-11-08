$ErrorActionPreference = "Stop"
$id = "1rxP6l5-BRZ1hjRzZpNH8r1SkH28Hk050"
$out = "activsg500_dataset.zip"
$base = "https://drive.google.com/uc?export=download&id=$id"
# first request: may be the file itself, or a virus-scan confirmation page
$resp = Invoke-WebRequest -Uri $base -UseBasicParsing -MaximumRedirection 10
if ($resp.Content -is [byte[]]) {
    [IO.File]::WriteAllBytes($out, $resp.Content)
    "downloaded $($resp.Content.Length) bytes"
} else {
    $html = $resp.Content
    # look for confirm token
    $m = [regex]::Match($html, 'name="confirm" value="([^"]+)"')
    $token = if ($m.Success) { $m.Groups[1].Value } else { $null }
    $m2 = [regex]::Match($html, 'confirm=([0-9A-Za-z_-]+)')
    if (-not $token -and $m2.Success) { $token = $m2.Groups[1].Value }
    if (-not $token) {
        # no confirm needed; maybe direct download link with uuid
        $m3 = [regex]::Match($html, 'href="(https://drive\.google\.com/uc[^\"]*export=download[^\"]*)"')
        if ($m3.Success) {
            $resp2 = Invoke-WebRequest -Uri $m3.Groups[1].Value -UseBasicParsing
            [IO.File]::WriteAllBytes($out, $resp2.Content)
            "downloaded $($resp2.Content.Length) bytes via link"
        } else {
            "NO DOWNLOAD LINK FOUND; page saved to activsg500_page.html"
            [IO.File]::WriteAllText("activsg500_page.html", $html)
        }
    } else {
        $url = "https://drive.google.com/uc?export=download&confirm=$token&id=$id"
        $resp2 = Invoke-WebRequest -Uri $url -UseBasicParsing
        [IO.File]::WriteAllBytes($out, $resp2.Content)
        "downloaded $($resp2.Content.Length) bytes with token"
    }
}
Get-Item $out | Select-Object Name, Length
