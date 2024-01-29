FROM node:18.15.0 as builder

RUN corepack enable

# Create app directory
WORKDIR /usr/src/app

# Install app dependencies
# A wildcard is used to ensure both package.json AND package-lock.json are copied
# where available (npm@5+)
COPY *.json ./
COPY yarn.lock ./
COPY .yarnrc.yml ./
COPY .yarn ./.yarn

# Bundle app source
COPY src/ ./src/

RUN yarn install && yarn run build && rm -rf node_modules && yarn workspaces focus --production

FROM node:18.15.0-alpine
COPY --from=builder /usr/src/app /usr/src/app

WORKDIR /usr/src/app
EXPOSE 5000
CMD [ "yarn", "run", "listen" ]
