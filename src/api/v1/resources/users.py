from fastapi import APIRouter, Depends, HTTPException, status, Security
from fastapi.security import OAuth2PasswordRequestForm

from src.api.v1.schemas import UserCreate, UserLogin, UserDetail, UserUpdate
from src.auth.auth import get_current_user, get_token
from src.services import UserService, get_user_service

router = APIRouter()


@router.post(path="/form-login",
             response_model=dict,
             summary="Авторизация пользователя через форму",
             tags=["auth"],

             )
def token(
        user_credentials: OAuth2PasswordRequestForm = Depends(),
        user_service: UserService = Depends(get_user_service),
) -> dict:
    """ Идентификация и аутентификация пользователя по данным из формы.
        Выдача пары access и refresh токенов
    """
    user = user_service.authenticate_user(user_credentials.username,
                                          user_credentials.password)
    if not user:
        # Если аутентификация не пройдена, отдаём 401 статус
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect login or password"
        )
    access_token = user_service.create_access_token(str(user.uuid))
    refresh_token = user_service.create_refresh_token(str(user.uuid))
    return {'access_token': access_token, 'refresh_token': refresh_token}


@router.post(
    path="/login",
    response_model=dict,
    summary="Авторизация пользователя",
    tags=["auth"],
)
def login(
        user_credentials: UserLogin,
        user_service: UserService = Depends(get_user_service),
) -> dict:
    """ Идентификация и аутентификация пользователя.
        Выдача пары access и refresh токенов
    """
    user = user_service.authenticate_user(user_credentials.username,
                                          user_credentials.password)
    if not user:
        # Если аутентификация не пройдена, отдаём 401 статус
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect login or password"
        )
    access_token = user_service.create_access_token(str(user.uuid))
    refresh_token = user_service.create_refresh_token(str(user.uuid))
    return {'access_token': access_token, 'refresh_token': refresh_token}


@router.post(
    path="/signup",
    response_model=UserDetail,
    summary="Создать пользователя",
    tags=["users"],
)
def user_create(
        user: UserCreate,
        user_service: UserService = Depends(get_user_service),
) -> UserDetail:
    """ Регистрация нового пользователя
    """
    if user_service.get_user_by_name(username=user.username):
        # Если пользователь с таким именем уже есть, отдаём 400 статус
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="User with this username already exist"
        )
    new_user: UserDetail = user_service.create_user(user=user)
    return new_user


@router.get(
    path="/users/me",
    response_model=UserDetail,
    summary="Информация о текущем пользователе пользователя",
    tags=["users"],
)
def user_get(
        user: str = Security(get_current_user),
        user_service: UserService = Depends(get_user_service),
) -> UserDetail:
    """ Возвращает данные о текущем пользователе
    """
    return user_service.get_user_detail(user)


@router.patch(
    path="/users/me",
    response_model=dict,
    summary="Изменить пользователя",
    tags=["users"],
)
def user_patch(
    user_changes: UserUpdate,
    user: str = Security(get_current_user),
    user_service: UserService = Depends(get_user_service),
) -> dict:
    updated_user = user_service.patch_user(user, user_changes)
    access_token = user_service.create_access_token(user)
    return {'user': updated_user, 'access_token': access_token}


@router.post(
    path="/logout",
    response_model=dict,
    summary="Выйти с одного устройства",
    tags=["auth"],
)
def logout(
        access_token: str = Security(get_token),
        user_service: UserService = Depends(get_user_service),
) -> dict:
    """ Выход пользователя. Добавление текущего access токена в черный список
    """
    user_service.add_token_to_black_list(access_token)
    return {"msg": "You have been logout"}


@router.post(
    path="/refresh",
    response_model=dict,
    summary="Обновление токенов",
    tags=["auth"],
)
def refresh(
        refresh_token: str = Security(get_token),
        user_service: UserService = Depends(get_user_service),
) -> dict:
    """ Обновление токенов
    """
    user = user_service.check_refresh_token(refresh_token)
    access_token = user_service.create_access_token(user)
    refresh_token = user_service.create_refresh_token(user)
    return {'access_token': access_token, 'refresh_token': refresh_token}


@router.post(
    path="/logout_all",
    response_model=dict,
    summary="Выйти со всех устройств",
    tags=["auth"],
)
def logout_all(
        access_token: str = Security(get_token),
        user_service: UserService = Depends(get_user_service),
) -> dict:
    """ Выход пользователя.
        Добавление текущего access токена в черный список,
        удаление всех refresh токенов пользователя
    """

    user = user_service.add_token_to_black_list(access_token)
    user_service.delete_all_refresh_tokens(user)
    return {"msg": "You have been logout"}
