from functools import lru_cache
from typing import Union
from fastapi import Depends, HTTPException, status
from sqlmodel import Session

from src.api.v1.schemas import UserCreate, UserDetail, UserUpdate
from src.auth.password import Password
from src.auth.token import Token
from src.core import config
from src.db import AbstractCache, get_cache, get_session
from src.db.cache import get_access_token_cache, get_refresh_token_cache
from src.models import User
from src.services import TokenCacheServiceMixin

__all__ = ("UserService", "get_user_service")


class UserService(TokenCacheServiceMixin):

    def get_user_detail(self, uuid: str) -> UserDetail:
        """
        Возвращает данные о пользователе(кроме пароля)
        """
        user: Union[User, None] = self.get_user_by_uuid(uuid)
        return UserDetail(**user.dict(exclude={'hashed_password'}))

    def get_user_by_name(self, username: str):
        """
        Выполняет поиск пользователя по username
        """
        return self.session.query(User).filter(User.username == username).first()

    def get_user_by_uuid(self, uuid: str):
        """
        Выполняет поиск пользователя по uuid
        """
        return self.session.query(User).filter(User.uuid == uuid).first()

    def authenticate_user(self, username: str, password: str):
        """
        Проверяет пару значений username, password на соответствие в БД
        """
        user = self.get_user_by_name(username)
        if not user:
            return False
        if not Password.verify_password(password, user.hashed_password):
            return False
        return user

    def create_access_token(self, user_uuid: str):
        """
        Для пользователя(uuid) создает access токен
        """
        access_token = Token.encode_token(user_uuid, 'access', config.ACCESS_TOKEN_EXPIRE_IN_SECONDS)
        return access_token

    def create_refresh_token(self, user_uuid: str):
        """
        Для пользователя(uuid) создает refresh токен
        """
        refresh_token = Token.encode_token(user_uuid, 'refresh', config.REFRESH_TOKEN_EXPIRE_IN_SECONDS)
        refresh_token_decode = Token.decode_token(refresh_token)
        self.refresh_token_cache.set(key=f"{refresh_token_decode['jti']}",
                                     value=f"{refresh_token_decode['sub']}",
                                     expire=config.REFRESH_TOKEN_EXPIRE_IN_SECONDS)

        return refresh_token

    def patch_user(self, uuid: str, patch_data: UserUpdate):
        """
        Изменяет данные пользователя
        """
        user: Union[User, None] = self.get_user_by_uuid(uuid)
        if not user:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                                detail="User doesn't exist",
                                headers={"WWW-Authenticate": "Bearer"}
                                )

        update_data = patch_data.dict(exclude_unset=True)
        if 'username' in update_data:
            if self.get_user_by_name(update_data['username']):
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail="User with this username already exist"
                )
        if 'password' in update_data:
            user.hashed_password = Password.get_password_hash(update_data['password'])
            del update_data['password']
        update_user = UserUpdate(**update_data)
        for var, value in vars(update_user).items():
            setattr(user, var, value) if value else None
        self.session.add(user)
        self.session.commit()
        self.session.refresh(user)
        return UserDetail(**user.dict(exclude={"hashed_password", }))

    def create_user(self, user: UserCreate) -> UserDetail:
        """Создает пользователя."""
        new_user = User(username=user.username, hashed_password=Password.get_password_hash(user.password),
                        email=user.email)
        self.session.add(new_user)
        self.session.commit()
        self.session.refresh(new_user)
        return UserDetail(**new_user.dict(exclude={"hashed_password", }))

    def add_token_to_black_list(self, access_token):
        """
        При выходе пользователя добавляет access токен в черный список
        """
        Token.check_token_type(access_token)
        if self.access_token_cache.get(key=f"{access_token['jti']}"):
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                                detail="Token is not valid. Please relogin",
                                headers={"WWW-Authenticate": "Bearer"}
                                )
        user = access_token['sub']
        if not user or not self.session.query(User.uuid).filter_by(uuid=user).first():
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                                detail="User doesn't exist",
                                headers={"WWW-Authenticate": "Bearer"}
                                )
        self.access_token_cache.set(key=f"{access_token['jti']}",
                                    value=f"{access_token['sub']}",
                                    expire=config.ACCESS_TOKEN_EXPIRE_IN_SECONDS)
        return access_token['sub']

    def check_refresh_token(self, refresh_token):
        """
        проверяет, что refresh токен валидный
        """
        Token.check_token_type(refresh_token, 'refresh')
        jti = refresh_token['jti']
        user = refresh_token['sub']
        if not user or not self.session.query(User.uuid).filter_by(uuid=user).first():
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                                detail="User doesn't exist",
                                headers={"WWW-Authenticate": "Bearer"}
                                )
        if self.refresh_token_cache.get(key=f"{jti}"):
            self.refresh_token_cache.delete(key=f"{jti}")
            return user
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                            detail="Token is not valid. Please relogin",
                            headers={"WWW-Authenticate": "Bearer"}
                            )

    def delete_all_refresh_tokens(self, user_uuid: str):
        """
        При выходе из всех устройств удаляем из кеша все refresh токены пользователя
        """
        for key in self.refresh_token_cache.scan_iter():
            if str(self.refresh_token_cache.get(key)) == user_uuid:
                self.refresh_token_cache.delete(key)
        return


# get_user_service — это провайдер UserService. Синглтон
@lru_cache()
def get_user_service(
        cache: AbstractCache = Depends(get_cache),
        session: Session = Depends(get_session),
        access_token_cache: AbstractCache = Depends(get_access_token_cache),
        refresh_token_cache: AbstractCache = Depends(get_refresh_token_cache),
) -> UserService:
    return UserService(cache=cache,
                       session=session,
                       access_token_cache=access_token_cache,
                       refresh_token_cache=refresh_token_cache)
