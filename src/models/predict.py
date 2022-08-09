from datetime import datetime

from sqlmodel import Field, SQLModel

__all__ = ("Predict",)


class Predict(SQLModel, table=True):
    id: int = Field(default=None, primary_key=True)
    ensemble_id: str = Field(nullable=False)
    exon_location: str = Field(nullable=False)
    created_at: datetime = Field(default=datetime.utcnow(), nullable=False)
