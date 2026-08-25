param(
  [ValidateSet('ALL','SUMMARY','Region_1','Region_2','Region_3','Region_4','Region_5','Region_6')]
  [string]$Mode = 'ALL',
  [string]$RunLabel = 'local_subset_validation'
)

$ErrorActionPreference = 'Stop'
$AnalysisRoot = 'D:\Xiaonan\CODEX_projects\Yanan_Xenium\colon_analysis'
$RepoRoot = Join-Path $AnalysisRoot 'YNH_Xenium_Colon'
$ProjectRoot = 'D:\Xiaonan\CODEX_projects\Yanan_Xenium'
$InputRoot = Join-Path $AnalysisRoot 'subset_input\colon_data'
$RunRoot = Join-Path $AnalysisRoot (Join-Path 'colon_qc_outputs' $RunLabel)
$TempRoot = Join-Path $AnalysisRoot (Join-Path 'tmp' $RunLabel)
$Rscript = 'D:\Programs\R-4.6.1\bin\Rscript.exe'
$RLibrary = 'D:\Programs\R_library'

foreach ($Path in @($AnalysisRoot, $RepoRoot, $InputRoot, $Rscript, $RLibrary)) {
  if (-not (Test-Path -LiteralPath $Path)) { throw "Required local path is missing: $Path" }
}
foreach ($WritablePath in @($RunRoot, $TempRoot)) {
  if (-not $WritablePath.StartsWith($AnalysisRoot, [System.StringComparison]::OrdinalIgnoreCase)) { throw "Writable path escaped colon_analysis: $WritablePath" }
  New-Item -ItemType Directory -Path $WritablePath -Force | Out-Null
}

$env:COLON_QC_MODE = 'LOCAL_SUBSET'
$env:COLON_PROJECT_ROOT = $ProjectRoot
$env:COLON_PIPELINE_REPO = $RepoRoot
$env:COLON_INPUT_ROOT = $InputRoot
$env:COLON_RUN_LABEL = $RunLabel
$env:COLON_RUN_ROOT = $RunRoot
$env:COLON_TEMP_ROOT = $TempRoot
$env:R_LIBS_USER = $RLibrary
$env:TMPDIR = $TempRoot
$env:TMP = $TempRoot
$env:TEMP = $TempRoot
$env:LC_ALL = 'C'

function Invoke-ColonNotebook([string]$NotebookName) {
  $Source = Join-Path (Join-Path $RepoRoot 'notebooks') $NotebookName
  $Executed = Join-Path (Join-Path $RunRoot 'executed_notebooks') $NotebookName
  & $Rscript (Join-Path $RepoRoot 'scripts\execute_r_notebook.R') $Source $Executed
  if ($LASTEXITCODE -ne 0) { throw "Notebook execution failed: $NotebookName" }
}

if ($Mode -eq 'ALL') {
  1..6 | ForEach-Object { Invoke-ColonNotebook ("01_QC_Region{0}.ipynb" -f $_) }
  Invoke-ColonNotebook '02_slide_QC_summary.ipynb'
} elseif ($Mode -eq 'SUMMARY') {
  Invoke-ColonNotebook '02_slide_QC_summary.ipynb'
} else {
  $Number = $Mode.Replace('Region_', '')
  Invoke-ColonNotebook ("01_QC_Region{0}.ipynb" -f $Number)
}
