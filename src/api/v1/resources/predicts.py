from http import HTTPStatus
from typing import Optional

from fastapi import APIRouter, Depends, HTTPException, Security

from src.api.v1.schemas import PredictCreate, PredictListResponse, PredictModel
from src.auth.auth import get_current_user
from src.services import PredictService, get_predict_service

router = APIRouter()


@router.get(
    path="/",
    response_model=PredictListResponse,
    summary="Список предсказаний",
    tags=["predicts"],
)
def predct_list(
    predict_service: PredictService = Depends(get_predict_service),
) -> PredictListResponse:
    predicts: dict = predict_service.get_predict_list()
    if not predicts:
        # Если предсказания не найдены, отдаём 404 статус
        raise HTTPException(status_code=HTTPStatus.NOT_FOUND,
                            detail="predicts are not found")
    return PredictListResponse(**predicts)


@router.get(
    path="/{predict_id}",
    response_model=PredictModel,
    summary="Получить определенное предсказание",
    tags=["predicts"],
)
def predict_detail(
    predict_id: int, predict_service: PredictService = Depends(get_predict_service),
) -> PredictModel:
    predict: Optional[dict] = predict_service.get_predict_detail(item_id=predict_id)
    if not predict:
        # Если предсказание не найдено, отдаём 404 статус
        raise HTTPException(status_code=HTTPStatus.NOT_FOUND,
                            detail="predict is not found")
    return PredictModel(**predict)


@router.post(
    path="/",
    response_model=PredictModel,
    summary="Создать предсказание",
    tags=["predicts"],
)
async def predict_create(
    predict: PredictCreate, user: str = Security(get_current_user),
    predict_service: PredictService = Depends(get_predict_service),
) -> PredictModel:
    predict: dict = await predict_service.create_predict(predict=predict)
    return PredictModel(**predict)
