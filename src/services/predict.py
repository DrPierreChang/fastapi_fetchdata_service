import json
from functools import lru_cache
from typing import Optional

from fastapi import Depends
from sqlmodel import Session

from src.api.v1.schemas import PredictCreate, PredictModel
from src.db import AbstractCache, get_cache, get_session
from src.models import Predict
from src.services import ServiceMixin

from ensemblrest import EnsemblRest

__all__ = ("PredictService", "get_predict_service")


async def fetch_transcript(transcript_id):
    ens_rest = EnsemblRest()
    transcript = ens_rest.getSequenceById(id=transcript_id,
                                                expand_5prime=5000,
                                                expand_3prime=5000)
    return f"{transcript}"


class PredictService(ServiceMixin):
    def get_predict_list(self) -> dict:
        """Получить список предсказаний."""
        predicts = self.session.query(Predict).order_by(Predict.created_at).all()
        return {"predicts": [PredictModel(**predict.dict())
                                for predict in predicts]}

    def get_predict_detail(self, item_id: int) -> Optional[dict]:
        """Получить детальную информацию предсказания."""
        if cached_predict := self.cache.get(key=f"{item_id}"):
            return json.loads(cached_predict)

        predict = self.session.query(Predict).filter(Predict.id == item_id).first()
        if predict:
            self.cache.set(key=f"{predict.id}", value=predict.json())
        return predict.dict() if predict else None

    async def create_predict(self, predict: PredictCreate) -> dict:
        """Создать предсказание."""
        exon_location = await fetch_transcript(predict.ensemble_id)
        new_predict = Predict(ensemble_id=predict.ensemble_id,
                              exon_location=exon_location)
        self.session.add(new_predict)
        self.session.commit()
        self.session.refresh(new_predict)
        return new_predict.dict()


# get_predict_service — это провайдер PredictService. Синглтон
@lru_cache()
def get_predict_service(
    cache: AbstractCache = Depends(get_cache),
    session: Session = Depends(get_session),
) -> PredictService:
    return PredictService(cache=cache, session=session)
