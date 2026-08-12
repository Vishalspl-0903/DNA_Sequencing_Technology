# Live view of the assembly batch.
#   Run in your own terminal:   powershell -File d:\DNA-Sequencing\work\watch.ps1
#   Ctrl+C to stop watching (this does NOT stop the batch).

$asm  = "d:\DNA-Sequencing\work\asm"
$task = "C:\Users\PREMN~1\AppData\Local\Temp\claude\d--DNA-Sequencing\24f8f323-e6fb-4a6c-abd9-a202e5e00ad6\tasks\bx2i25nen.output"
$total = 49   # 48 cohort cells + the sim01_d50 pilot

while ($true) {
    Clear-Host
    $done = (Get-ChildItem $asm -Directory -ErrorAction SilentlyContinue |
             Where-Object { Test-Path "$($_.FullName)\assembly_graph.gfa" }).Count
    $left = (wsl -d Ubuntu -- bash -c "ls -d /home/premn/scratch/*/ 2>/dev/null | wc -l").Trim()

    $pct = [math]::Round(100 * $done / $total)
    $barLen = 40
    $fill = [math]::Round($barLen * $done / $total)
    $bar = ("#" * $fill).PadRight($barLen, ".")

    Write-Host ""
    Write-Host "  PLASMID PIPELINE - assembly batch" -ForegroundColor Cyan
    Write-Host "  ---------------------------------"
    Write-Host "  [$bar] $pct%"
    Write-Host "  graphs built : $done / $total"
    Write-Host "  queued       : $left"
    Write-Host "  time         : $(Get-Date -Format 'HH:mm:ss')"
    Write-Host ""

    # per-arm timing from the manifests written so far
    $man = Get-ChildItem $asm -Directory -ErrorAction SilentlyContinue |
        ForEach-Object { $p = "$($_.FullName)\manifest.json"
                         if (Test-Path $p) { Get-Content $p -Raw | ConvertFrom-Json } }
    if ($man) {
        Write-Host "  mean Flye seconds by arm:" -ForegroundColor DarkGray
        $man | Group-Object arm | ForEach-Object {
            $avg = [math]::Round(($_.Group | Measure-Object seconds_flye -Average).Average)
            Write-Host ("    {0,-8} n={1,-3} avg={2}s" -f $_.Name, $_.Count, $avg)
        }
        Write-Host ""
    }

    Write-Host "  recent:" -ForegroundColor DarkGray
    if (Test-Path $task) {
        Get-Content $task -Tail 12 | ForEach-Object { Write-Host "    $_" }
    }
    Write-Host ""
    Write-Host "  (Ctrl+C stops watching, not the batch)" -ForegroundColor DarkGray
    Start-Sleep -Seconds 10
}
