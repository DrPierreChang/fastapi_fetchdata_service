from sqlmodel import Session

from src.db import AbstractCache


class ServiceMixin:
    def __init__(self, cache, session):
        self.cache: AbstractCache = cache
        self.session: Session = session


class TokenCacheServiceMixin(ServiceMixin):
    def __init__(self, cache: AbstractCache,
                 session: Session,
                 access_token_cache: AbstractCache,
                 refresh_token_cache: AbstractCache):
        super().__init__(cache, session)
        self.access_token_cache: AbstractCache = access_token_cache
        self.refresh_token_cache: AbstractCache = refresh_token_cache
