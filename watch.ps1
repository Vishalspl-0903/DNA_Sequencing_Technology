# Live view of the cohort assembly batch.
#
#   powershell -File D:\Plasmid-GNN\DNA_Sequencing_Technology\work\watch.ps1
#
# Ctrl+C stops WATCHING. It does not stop the batch -- the batch runs inside
# WSL and is not a child of this window.
#
# Rewritten 2026-08-14 for cohort v3. The previous version pointed at
# d:\DNA-Sequencing (the tree moved), hardcoded $total = 49, and tailed a task
# output file under a stale session UUID belonging to a different Windows user.
# It now derives everything from disk, so it cannot go stale that way again.

$root  = "D:\Plasmid-GNN\DNA_Sequencing_Technology\work"
$asm   = Join-Path $root "asm"
$log   = Join-Path $root "logs\cohort_v3.log"
$sim   = Join-Path $root "sim"
$every = 30      # seconds between refreshes

# total cells = one reference FASTA per isolate x {lig,rap} x {30x,15x}
$isolates = (Get-ChildItem $sim -Filter "sim*_ref.fasta" -ErrorAction SilentlyContinue).Count
if ($isolates -eq 0) { Write-Host "No references in $sim -- run 02_build_sim_refs.py first."; exit 1 }
$total = $isolates * 4

$start = Get-Date
while ($true) {
    $done = 0
    if (Test-Path $asm) {
        $done = (Get-ChildItem $asm -Directory -ErrorAction SilentlyContinue |
                 Where-Object { (Test-Path (Join-Path $_.FullName "assembly_graph.gfa")) -and
                                ((Get-Item (Join-Path $_.FullName "assembly_graph.gfa")).Length -gt 0) }).Count
    }
    $left = 0
    try {
        $r = wsl -d Ubuntu -- bash -c "ls -d ~/scratch/*/ 2>/dev/null | wc -l"
        $left = [int]($r | Select-Object -First 1).ToString().Trim()
    } catch { }

    $elapsed = (Get-Date) - $start
    $pct  = if ($total) { [math]::Round(100 * $done / $total, 1) } else { 0 }
    $fill = if ($total) { [math]::Round(40 * $done / $total) } else { 0 }
    $bar  = ("#" * $fill).PadRight(40, ".")

    Clear-Host
    Write-Host ""
    Write-Host "  COHORT v3 -- $isolates isolates x 4 cells" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "  [$bar] $done / $total  ($pct%)"
    Write-Host ""
    Write-Host ("  read sets staged in scratch : {0}" -f $left)
    Write-Host ("  watching for               : {0:hh\:mm\:ss}" -f $elapsed)

    # rate and ETA come from THIS window's observations, so they are only
    # meaningful once a few assemblies have landed since it started
    if ($done -gt 0 -and $elapsed.TotalMinutes -gt 2) {
        $rate = $done / $elapsed.TotalMinutes
        if ($rate -gt 0) {
            $etaMin = ($total - $done) / $rate
            Write-Host ("  rate                       : {0:N2} assemblies/min" -f $rate)
            Write-Host ("  ETA (this window's rate)   : {0:N0} min  (~{1:N1} h)" -f $etaMin, ($etaMin / 60))
        }
    } else {
        Write-Host "  rate                       : (need ~2 min of observation)"
    }

    if (Test-Path $log) {
        Write-Host ""
        Write-Host "  --- last 12 log lines ---" -ForegroundColor DarkGray
        Get-Content $log -Tail 12 -ErrorAction SilentlyContinue |
            ForEach-Object { Write-Host ("  " + $_.TrimEnd()) }
    } else {
        Write-Host ""
        Write-Host "  (no $log yet -- batch not started)" -ForegroundColor DarkGray
    }

    Write-Host ""
    Write-Host "  refreshing every $every s -- Ctrl+C stops watching, not the batch" -ForegroundColor DarkGray
    Start-Sleep -Seconds $every
}
