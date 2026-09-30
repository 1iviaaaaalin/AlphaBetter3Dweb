from pathlib import Path
from tempfile import TemporaryDirectory

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles

from utils import read_residue_scores, run_calc_stability, write_b_factor

BASE_DIR = Path(__file__).resolve().parent
FRONTEND_DIR = BASE_DIR / "frontend" / "dist"
SCRIPT_PATH = BASE_DIR / "calc_stability.py"
DI_DB = BASE_DIR / "data" / "db_dipeptides.csv"
PENTA_DB = BASE_DIR / "data" / "db_pentapeptides.csv"

app = FastAPI(title="AlphaBetter", version="1.0.0")


@app.get("/api/health")
def health():
    return {"status": "ok", "service": "AlphaBetter"}


@app.post("/api/upload")
async def upload_pdb(file: UploadFile = File(...)):
    filename = Path(file.filename or "protein.pdb").name
    if not filename.lower().endswith(".pdb"):
        raise HTTPException(status_code=400, detail="只支持 .pdb 文件")

    content = await file.read()
    if not content:
        raise HTTPException(status_code=400, detail="上传的 PDB 文件为空")
    if len(content) > 20 * 1024 * 1024:
        raise HTTPException(status_code=413, detail="文件过大，请上传不超过 20 MB 的 PDB 文件")

    with TemporaryDirectory(prefix="alphabetteR_") as temp_dir:
        temp = Path(temp_dir)
        raw_pdb = temp / filename
        score_dir = temp / "scores"
        output_pdb = temp / f"colored_{filename}"
        raw_pdb.write_bytes(content)

        try:
            score_csv = run_calc_stability(SCRIPT_PATH, raw_pdb, DI_DB, PENTA_DB, score_dir)
            score_map = read_residue_scores(score_csv)
            write_b_factor(raw_pdb, output_pdb, score_map)
        except Exception as exc:
            raise HTTPException(status_code=500, detail=str(exc)) from exc

        result = output_pdb.read_bytes()

    download_name = f"colored_{filename}"
    return Response(
        content=result,
        media_type="chemical/x-pdb",
        headers={
            "Content-Disposition": f'attachment; filename="{download_name}"',
            "X-AlphaBetter-Status": "ok",
        },
    )


# 单一网址部署：FastAPI 同时负责 API 和 React 静态页面。
if FRONTEND_DIR.exists():
    app.mount("/", StaticFiles(directory=FRONTEND_DIR, html=True), name="frontend")
else:
    @app.get("/")
    def missing_frontend():
        return {"message": "AlphaBetter 后端已启动，但前端尚未构建。"}
