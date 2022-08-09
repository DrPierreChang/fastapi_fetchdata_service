from fastapi import Security, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordBearer
from sqlmodel import Session

from src.auth.token import Token
from src.db import get_session, get_access_token_cache, AbstractCache
from src.models import User

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/api/v1/form-login")


def get_token(token: str = Security(oauth2_scheme)):
    """ Получение токена
    """
    payload = Token.decode_token(token)
    return payload


def get_current_user(payload: dict = Security(get_token),
                     session: Session = Depends(get_session),
                     cache: AbstractCache = Depends(get_access_token_cache)):
    """ Проверка uuid пользователя
    """
    Token.check_token_type(payload)
    if cache.get(key=f"{payload['jti']}"):
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                            detail="Token is not valid. Please relogin",
                            headers={"WWW-Authenticate": "Bearer"}
                            )
    user = payload['sub']
    if not user or not session.query(User.uuid).filter_by(uuid=user).first():
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                            detail="User doesn't exist",
                            headers={"WWW-Authenticate": "Bearer"}
                            )
    return user
