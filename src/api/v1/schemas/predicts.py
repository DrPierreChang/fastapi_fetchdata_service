from datetime import datetime

from pydantic import BaseModel

__all__ = (
    "PredictModel",
    "PredictCreate",
    "PredictListResponse",
)


class PredictBase(BaseModel):
    ensemble_id: str
    exon_location: str


class PredictCreate(PredictBase):
    ...


class PredictModel(PredictBase):
    id: int
    created_at: datetime


class PredictListResponse(BaseModel):
    predicts: list[PredictModel] = []
